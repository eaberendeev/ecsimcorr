#include "simulation_symplectic.h"

#include <cmath>
#include <exception>
#include <filesystem>
#include <stdexcept>

#include "config/boundary_config.h"
#include "config/domain_config.h"
#include "config/symplectic_config.h"
#include "diagnostics/DiagnosticMeasurements.h"
#include "diagnostics/OutputFactory.h"
#include "io/recovery.h"
#include "runtime/profiler.h"

void SimulationSymplectic::init() {
    RECORD_TIMER;
    validate_symplectic_config(system_config, particles_config);

    // Check the actual handler semantics, including omitted periodic axes.
    Domain checked_domain = domain_from_config(system_config);
    BoundaryConditionHandler checked_bc = boundary_conditions_from_config(system_config, checked_domain);
    for (int a = 0; a < 3; ++a)
        if (!checked_bc.is_periodic(a))
            throw std::invalid_argument("symplectic requires three periodic axes");

    if (system_config.value("StartFromTime", 0) > 0) {
        for (const char* path : {"../Recovery/Fields/FieldE.backup", "../Recovery/Fields/FieldB.backup"})
            if (!std::filesystem::is_regular_file(path))
                throw std::runtime_error("missing symplectic recovery field");
        for (const nlohmann::json& sp : particles_config.at("particles")) {
            const std::string name = sp.at("Name");
            if (!std::filesystem::is_regular_file("../Recovery/Particles/" + name + "/" + name + ".backup"))
                throw std::runtime_error("missing symplectic recovery species: " + name);
        }
    }

    Simulation::init();
    validate_initial_state();
    for (const Species::value_type& entry : species) entry.second->currentOnGrid.setZero();

    Field3d rho(domain.size(), 1);
    collect_charge_density(rho);
    gauss_initial_.resize(domain.size(), 1);
    const int g = GHOST_CELLS;
    const Vector3I n = domain.num_cells();
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k) {
                if (std::abs(symplectic::divergence(domain, magnetic, i, j, k, true)) > 1e-12)
                    throw std::invalid_argument("symplectic requires divergence-free initial B");
                gauss_initial_(i + g, j + g, k + g, 0) =
                    symplectic::divergence(domain, electric, i, j, k, false) - rho(i + g, j + g, k + g, 0);
            }

    energy_initial_ = symplectic::field_energy(domain, electric) + symplectic::field_energy(domain, magnetic);
    for (const Species::value_type& entry : species) energy_initial_ += entry.second->get_kinetic_energy();
}

void SimulationSymplectic::validate_initial_state() const {
    RECORD_TIMER;
    if (!electric.data().allFinite() || !magnetic.data().allFinite())
        throw std::runtime_error("non-finite symplectic fields");
    for (const Species::value_type& entry : species)
        for (int cell = 0; cell < entry.second->size(); ++cell)
            for (const Particle& p : (*entry.second)(cell)) {
                for (int d = 0; d < 3; ++d)
                    if (!std::isfinite(p.coord[d]) || !std::isfinite(p.velocity[d]))
                        throw std::runtime_error("non-finite symplectic particle");
                if (!domain.contains(p.coord))
                    throw std::runtime_error("symplectic particle outside periodic box");
            }
}

void SimulationSymplectic::init_fields() {
    RECORD_TIMER;
    electric = Field3d::Zero(domain.size(), 3);
    magnetic = Field3d::Zero(domain.size(), 3);
    zero_field_ = Field3d::Zero(domain.size(), 3);
    if (system_config.value("StartFromTime", 0) > 0) {
        read_fields_from_recovery(electric, magnetic);
    } else {
        const nlohmann::json cfg = system_config.value("Symplectic", nlohmann::json::object());
        const Vector3R e = initial_symplectic_vector(cfg, "initial_electric"),
                       b = initial_symplectic_vector(cfg, "initial_magnetic");
#pragma omp parallel for schedule(static)
        for (int i = 0; i < electric.capacity(); ++i) {
            electric[i] = e[i % 3];
            magnetic[i] = b[i % 3];
        }
    }
    symplectic::fill_periodic_ghosts(domain, electric);
    symplectic::fill_periodic_ghosts(domain, magnetic);
}

int SimulationSymplectic::reset_deposit_buffers(std::vector<Field3d>& buffers, int components,
                                                const ParticlesArray& sp) const {
    RECORD_TIMER;
    const int threads = sp.get_total_num_of_particles() >= 4096 ? omp_get_max_threads() : 1;
    if (static_cast<int>(buffers.size()) < threads)
        buffers.resize(threads);
#pragma omp parallel for num_threads(threads) if (threads > 1) schedule(static)
    for (int t = 0; t < threads; ++t) {
        if (buffers[t].sizes() != domain.size() || buffers[t].nd() != components)
            buffers[t].resize(domain.size(), components);
        buffers[t].setZero();
    }
    return threads;
}

void SimulationSymplectic::deposit_charge_density(const ParticlesArray& sp, Field3d& rho) const {
    RECORD_TIMER;
    const double weight = sp.is_neutral() ? sp.mpw() : sp.charge() * sp.mpw();
    const int threads = reset_deposit_buffers(density_buffers_, 1, sp);
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        Field3d& local = density_buffers_[omp_get_thread_num()];
#pragma omp for schedule(static)
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const Particle& p : sp(cell)) symplectic::deposit_charge(domain, local, p.coord, weight);
    }
    if (rho.sizes() != domain.size() || rho.nd() != 1)
        rho.resize(domain.size(), 1);
    rho.setZero();
#pragma omp parallel for num_threads(threads) if (threads > 1) schedule(static)
    for (int i = 0; i < rho.capacity(); ++i) {
        double value = 0;
        for (int t = 0; t < threads; ++t) value += density_buffers_[t][i];
        rho[i] = value;
    }
    symplectic::fill_periodic_ghosts(domain, rho);
}

void SimulationSymplectic::collect_charge_density(Field3d& field) const {
    RECORD_TIMER;
    field.resize(domain.size(), 1);
    field.setZero();
    for (const Species::value_type& entry : species) {
        if (entry.second->is_neutral())
            continue;
        Field3d rho(domain.size(), 1);
        deposit_charge_density(*entry.second, rho);
        field += rho;
    }
}

void SimulationSymplectic::prepare_step(int) {
    RECORD_TIMER;
    for (const Species::value_type& entry : species) entry.second->begin_step();
}

void SimulationSymplectic::kinetic_flow(bool reverse, double dt) {
    RECORD_TIMER;
    for (const Species::value_type& entry : species) {
        const bool isNeutral = entry.second->is_neutral();
        const double qm = isNeutral ? 0 : entry.second->charge() / entry.second->mass();
        const double qw = isNeutral ? 0 : entry.second->charge() * entry.second->mpw();
        int threads;
        if (isNeutral) {
            threads = entry.second->get_total_num_of_particles() >= 4096 ? omp_get_max_threads() : 1;
        } else {
            timer::commonTimer buffers_timer("symplectic current buffers");
            threads = reset_deposit_buffers(current_buffers_, 3, *entry.second);
        }
        std::exception_ptr failure;
        timer::commonTimer particle_timer("symplectic particle flow");
#pragma omp parallel num_threads(threads) if (threads > 1)
        {
            Field3d* integral = isNeutral ? nullptr : &current_buffers_[omp_get_thread_num()];
#pragma omp for schedule(static)
            for (int cell = 0; cell < entry.second->size(); ++cell) {
                try {
                    if (isNeutral) {
                        for (Particle& particle : (*entry.second)(cell))
                            symplectic::drift_neutral(domain, particle, dt);
                    } else {
                        symplectic::advance_bucket(domain, (*entry.second)(cell), cell, reverse, dt, qm, qw, magnetic,
                                                   *integral);
                    }
                } catch (...) {
#pragma omp critical(symplectic_trajectory_failure)
                    {
                        if (!failure)
                            failure = std::current_exception();
                    }
                }
            }
        }
        particle_timer.finish();
        if (failure)
            std::rethrow_exception(failure);
        if (isNeutral)
            continue;
        // H_axis holds B fixed; E does not enter this particle flow. Deposits
        // therefore commute and can be summed before updating the common E.
        timer::commonTimer reduction_timer("symplectic current reduction");
#pragma omp parallel for num_threads(threads) if (threads > 1) schedule(static)
        for (int i = 0; i < electric.capacity(); ++i) {
            double integral = 0;
            for (int t = 0; t < threads; ++t) integral += current_buffers_[t][i];
            electric[i] -= integral;
            entry.second->currentOnGrid[i] += integral;
        }
        reduction_timer.finish();
    }
}

void SimulationSymplectic::electric_flow(double dt) {
    RECORD_TIMER;
    symplectic::electric_flow(domain, magnetic, electric, dt);
    for (const Species::value_type& entry : species) {
        if (entry.second->is_neutral())
            continue;
        const double impulse = dt * entry.second->charge() / entry.second->mass();
        const int threads = entry.second->get_total_num_of_particles() >= 4096 ? omp_get_max_threads() : 1;
#pragma omp parallel for num_threads(threads) if (threads > 1) schedule(static)
        for (int cell = 0; cell < entry.second->size(); ++cell)
            for (Particle& p : (*entry.second)(cell))
                p.velocity += impulse * symplectic::gather<false>(domain, electric, p.coord);
    }
}

void SimulationSymplectic::make_step(int) {
    RECORD_TIMER;
    // Xiao et al. (2015), Eq. (55): palindrome of exact Hamiltonian flows.
    // Unlike ECSIM, no half-time particle positions or averaged E are stored.
    kinetic_flow(false, 0.5 * dt_);
    symplectic::magnetic_flow(domain, electric, magnetic, 0.5 * dt_);
    electric_flow(dt_);
    symplectic::magnetic_flow(domain, electric, magnetic, 0.5 * dt_);
    kinetic_flow(true, 0.5 * dt_);
    symplectic::fill_periodic_ghosts(domain, electric);
    symplectic::fill_periodic_ghosts(domain, magnetic);
    for (const Species::value_type& entry : species) {
        if (!entry.second->is_neutral()) {
            entry.second->currentOnGrid *= 1.0 / dt_;
            symplectic::fill_periodic_ghosts(domain, entry.second->currentOnGrid);
        }
        // Each block scans source cells in order. Insert only after every block
        // finishes, so a particle cannot be advanced twice or race a target bucket.
        timer::commonTimer rebin_timer("symplectic rebin");
        constexpr int cells_per_block = 512;
        const int blocks = 1 + (entry.second->size() - 1) / cells_per_block;
        moved_buffers_.resize(blocks);
#pragma omp parallel for schedule(static)
        for (int block = 0; block < blocks; ++block) {
            std::vector<Particle>& moved = moved_buffers_[block];
            moved.clear();
            const int begin = block * cells_per_block;
            const int end = begin + std::min(cells_per_block, entry.second->size() - begin);
            for (int cell = begin; cell < end; ++cell) {
                std::vector<Particle>& bucket = (*entry.second)(cell);
                size_t p = 0;
                while (p < bucket.size()) {
                    const Vector3I target = domain.storage_cell_index(bucket[p].coord);
                    if (domain.sind(target[0], target[1], target[2]) == cell) {
                        ++p;
                    } else {
                        moved.push_back(bucket[p]);
                        bucket[p] = bucket.back();
                        bucket.pop_back();
                    }
                }
            }
        }
        for (std::vector<Particle>& moved : moved_buffers_) {
            entry.second->add_particles(moved);
            moved.clear();
        }
        rebin_timer.finish();
    }
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
    double kinetic = 0;
    for (const Species::value_type& entry : species) {
        const double energy = entry.second->get_kinetic_energy();
        diagnostic_->addEnergy(entry.first, energy);
        kinetic += energy;
    }
    const double e = symplectic::field_energy(domain, electric), b = symplectic::field_energy(domain, magnetic);
    diagnostic_->addEnergy("energyFieldE", e);
    diagnostic_->addEnergy("energyFieldB", b);
    if (!std::isfinite(kinetic + e + b))
        throw std::runtime_error("non-finite symplectic energy");
    diagnostic_->addEnergy("energyTotal", kinetic + e + b);
    diagnostic_->addEnergy("energyError", kinetic + e + b - energy_initial_);
    Field3d rho(domain.size(), 1);
    collect_charge_density(rho);
    double gauss_drift = 0, gauss_error = 0, div_b = 0;
    const Vector3I n = domain.num_cells();
    const int g = GHOST_CELLS;
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k) {
                const double defect =
                    symplectic::divergence(domain, electric, i, j, k, false) - rho(i + g, j + g, k + g, 0);
                gauss_error = std::max(gauss_error, std::abs(defect));
                gauss_drift = std::max(gauss_drift, std::abs(defect - gauss_initial_(i + g, j + g, k + g, 0)));
                div_b = std::max(div_b, std::abs(symplectic::divergence(domain, magnetic, i, j, k, true)));
            }
    diagnostic_->addEnergy("gaussError", gauss_error);
    diagnostic_->addEnergy("gaussDrift", gauss_drift);
    diagnostic_->addEnergy("divB", div_b);
    write_diagnostics(timestep, *diagnostic_,
                      {domain, species, electric, magnetic, zero_field_, zero_field_, bc_handler, SHAPE}, outputs_);
}

void SimulationSymplectic::finalize() {
    RECORD_TIMER;
    for (std::unique_ptr<IDiagnosticOutput>& out : outputs_) out->finalize();
    Simulation::finalize();
}
