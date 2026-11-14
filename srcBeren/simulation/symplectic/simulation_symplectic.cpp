#include "simulation_symplectic.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "config/symplectic_config.h"
#include "diagnostics/DensityMeasurements.h"
#include "diagnostics/DiagnosticMeasurements.h"
#include "diagnostics/OutputFactory.h"
#include "io/recovery.h"
#include "runtime/profiler.h"
#include "whitney_pic.h"

namespace {
// The kernels index by position without checks: reject bad input once, before stepping.
void validate_initial_state(const Domain& domain, const Species& species, const Field3d& electric,
                            const Field3d& magnetic) {
    RECORD_TIMER;
    if (!electric.data().allFinite() || !magnetic.data().allFinite())
        throw std::runtime_error("non-finite symplectic fields");

    for (const Species::value_type& entry : species) {
        const ParticlesArray& sp = *entry.second;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const Particle& p : sp(cell)) {
                for (int d = 0; d < 3; ++d)
                    if (!std::isfinite(p.coord[d]) || !std::isfinite(p.velocity[d]))
                        throw std::runtime_error("non-finite symplectic particle");
                if (!domain.contains(p.coord))
                    throw std::runtime_error("symplectic particle outside periodic box");
            }
    }
}
}   // namespace

void SimulationSymplectic::init() {
    RECORD_TIMER;
    validate_symplectic_config(system_config, particles_config);
    Simulation::init();
    validate_initial_state(domain, species, electric, magnetic);

    if (symplectic::max_magnetic_divergence(domain, magnetic) > 1e-12)
        throw std::invalid_argument("symplectic requires divergence-free initial B");

    for (const Species::value_type& entry : species) entry.second->currentOnGrid.setZero();
    gauss_initial_ = gauss_defect();
    energy_initial_ = field_energy(electric) + field_energy(magnetic);
    for (const Species::value_type& entry : species) energy_initial_ += entry.second->get_kinetic_energy();
}

void SimulationSymplectic::init_fields() {
    RECORD_TIMER;
    // The handler also knows the axes that the configuration omits.
    for (int axis = 0; axis < 3; ++axis)
        if (!bc_handler.is_periodic(axis))
            throw std::invalid_argument("symplectic requires three periodic axes");

    // The periodic Yee operators update every stored image of a node alike,
    // so the fields never need a separate ghost exchange.
    maxwell.init(assemble_maxwell_stencils(domain, bc_handler));

    electric = Field3d::Zero(domain.size(), 3);
    magnetic = Field3d::Zero(domain.size(), 3);
    zero_field_ = Field3d::Zero(domain.size(), 3);
    scratch_ = Field3d::Zero(domain.size(), 3);

    if (system_config.value("StartFromTime", 0) > 0) {
        read_fields_from_recovery(electric, magnetic);
        return;
    }
    const nlohmann::json config = system_config.value("Symplectic", nlohmann::json::object());
    const Vector3R e = initial_symplectic_vector(config, "initial_electric");
    const Vector3R b = initial_symplectic_vector(config, "initial_magnetic");
#pragma omp parallel for schedule(static)
    for (int i = 0; i < electric.capacity(); ++i) {
        electric[i] = e[i % 3];
        magnetic[i] = b[i % 3];
    }
}

double SimulationSymplectic::field_energy(const Field3d& field) const {
    return 0.5 * dot_product_sum(field, field, bc_handler.active_range(domain.grid));
}

Field3d SimulationSymplectic::gauss_defect() const {
    RECORD_TIMER;
    Field3d defect(domain.size(), 1);
    spmv(maxwell.divE, electric, defect);
    // The Whitney nodal charge is the Linear shape of the shared density deposit.
    defect -= measure_charge_density(species, domain, bc_handler, ShapeType::Linear);
    return defect;
}

void SimulationSymplectic::prepare_step(int) {
    RECORD_TIMER;
    for (const Species::value_type& entry : species) entry.second->currentOnGrid.setZero();
}

void SimulationSymplectic::kinetic_flow(bool reverse, double dt) {
    RECORD_TIMER;
    for (const Species::value_type& entry : species) {
        ParticlesArray& sp = *entry.second;
        if (sp.is_neutral()) {
            sp.move(dt);
        } else {
            Field3d& current = scratch_;
            current.setZero();
            if (!symplectic::kinetic_flow(sp, magnetic, dt, reverse, current)) [[unlikely]] {
                std::cerr << "symplectic: " << sp.name() << " moves more than a cell per axis; reduce Dt" << std::endl;
                std::abort();
            }
            bc_handler.apply_to_fields(current, FieldType::CURRENT, domain);

            // E does not enter the kinetic flows, so the deposits of all particles commute.
            electric -= current;
            blas::axpby(1.0 / dt_, current, 1.0, sp.currentOnGrid);
        }
        // The next flow needs every particle inside the box and in the bucket of its cell.
        bc_handler.apply_to_particles(sp, species, domain);
        sp.update_cells(domain);
    }
}

void SimulationSymplectic::electric_flow(double dt) {
    RECORD_TIMER;
    spmv(maxwell.curlE, electric, scratch_);
    blas::axpby(-dt, scratch_, 1.0, magnetic);
    for (const Species::value_type& entry : species) {
        ParticlesArray& sp = *entry.second;
        if (!sp.is_neutral())
            symplectic::kick_velocities(sp, electric, dt * sp.charge() / sp.mass());
    }
}

void SimulationSymplectic::magnetic_flow(double dt) {
    RECORD_TIMER;
    spmv(maxwell.curlB, magnetic, scratch_);
    blas::axpby(dt, scratch_, 1.0, electric);
}

void SimulationSymplectic::make_step(int) {
    RECORD_TIMER;
    // Xiao et al. (2015), Eq. (55): palindrome of exact Hamiltonian flows.
    // Unlike ECSIM, no half-time particle positions or averaged E are stored.
    kinetic_flow(false, 0.5 * dt_);
    magnetic_flow(0.5 * dt_);
    electric_flow(dt_);
    magnetic_flow(0.5 * dt_);
    kinetic_flow(true, 0.5 * dt_);
}

void SimulationSymplectic::make_diagnostic(int timestep) {
    RECORD_TIMER;
    if (timestep == 0)
        timestep = system_config.value("StartTimeStep", 0);
    if (!diagnostic_) {
        diagnostic_ = std::make_unique<Diagnostics>(species);
        outputs_ = OutputFactory::create(system_config.value("diagnostics", nlohmann::json::object()), domain, species,
                                         system_config, electric, magnetic, zero_field_, zero_field_);
    }

    double energy = 0;
    for (const Species::value_type& entry : species) {
        const double kinetic = entry.second->get_kinetic_energy();
        diagnostic_->addEnergy(entry.first, kinetic);
        diagnostic_->addEnergy(entry.first + "Particles",
                               static_cast<double>(entry.second->get_total_num_of_particles()));
        energy += kinetic;
    }
    const double e = field_energy(electric), b = field_energy(magnetic);
    energy += e + b;
    diagnostic_->addEnergy("energyFieldE", e);
    diagnostic_->addEnergy("energyFieldB", b);
    diagnostic_->addEnergy("energyTotal", energy);
    diagnostic_->addEnergy("energyError", energy - energy_initial_);

    const Field3d defect = gauss_defect();
    diagnostic_->addEnergy("gaussError", defect.data().lpNorm<Eigen::Infinity>());
    diagnostic_->addEnergy("gaussDrift", (defect - gauss_initial_).data().lpNorm<Eigen::Infinity>());
    diagnostic_->addEnergy("divB", symplectic::max_magnetic_divergence(domain, magnetic));

    write_diagnostics(timestep, *diagnostic_,
                      {domain, species, electric, magnetic, zero_field_, zero_field_, bc_handler, SHAPE}, outputs_);
}

void SimulationSymplectic::finalize() {
    RECORD_TIMER;
    for (std::unique_ptr<IDiagnosticOutput>& out : outputs_) out->finalize();
    Simulation::finalize();
}
