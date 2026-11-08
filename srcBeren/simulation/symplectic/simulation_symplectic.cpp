#include "simulation_symplectic.h"

#include <cmath>
#include <exception>
#include <filesystem>
#include <stdexcept>

#include "config/boundary_config.h"
#include "config/domain_config.h"
#include "diagnostics/DiagnosticMeasurements.h"
#include "diagnostics/OutputFactory.h"
#include "io/recovery.h"
#include "runtime/profiler.h"

namespace {
Vector3R initial_vector(const nlohmann::json& config, const char* name) {
    const auto values = config.value(name, std::vector<double>{0, 0, 0});
    if (values.size() != 3)
        throw std::invalid_argument(std::string(name) + " requires three components");
    for (double x : values)
        if (!std::isfinite(x))
            throw std::invalid_argument(std::string(name) + " must be finite");
    return {values[0], values[1], values[2]};
}
}   // namespace

void SimulationSymplectic::init() {
    RECORD_TIMER;
    // Validate the Hamiltonian model before allocating the mesh or loading particles.
    const double dt = system_config.at("Dt").get<double>();
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("symplectic requires finite positive Dt");
    double inverse_h2 = 0;
    for (const char* key : {"Dx", "Dy", "Dz"}) {
        const double h = system_config.at(key).get<double>();
        if (!std::isfinite(h) || h <= 0)
            throw std::invalid_argument("symplectic requires positive cell sizes");
        inverse_h2 += 1 / (h * h);
    }
    if (dt * std::sqrt(inverse_h2) >= 1)
        throw std::invalid_argument("symplectic explicit Maxwell CFL requires Dt*sqrt(sum(1/h^2)) < 1");
    for (const char* key : {"NumCellsX", "NumCellsY", "NumCellsZ"})
        if (system_config.at(key).get<int>() < 1)
            throw std::invalid_argument("symplectic requires nonempty grid");
    if (system_config.contains("CylinderDomain") || system_config.value("DampingType", "None") != "None" ||
        system_config.value("Collider", "None") != "None" ||
        (system_config.contains("Collisions") && !system_config.at("Collisions").empty()))
        throw std::invalid_argument("symplectic supports collisionless periodic boxes without damping");
    for (const char* key : {"ExternalFieldE", "ExternalFieldB"})
        if (system_config.contains(key) && !system_config.at(key).empty())
            throw std::invalid_argument(
                "symplectic: use Symplectic.initial_electric/initial_magnetic, not external fields");
    for (const auto& boundary : system_config.at("Boundary_conditions"))
        if (boundary.size() != 1 || !boundary.contains("periodic"))
            throw std::invalid_argument("symplectic supports only periodic boundary conditions");
    const auto config = system_config.value("Symplectic", nlohmann::json::object());
    for (auto it = config.begin(); it != config.end(); ++it)
        if (it.key() != "initial_electric" && it.key() != "initial_magnetic")
            throw std::invalid_argument("unknown Symplectic option: " + it.key());
    initial_vector(config, "initial_electric");
    initial_vector(config, "initial_magnetic");
    for (const auto& sp : particles_config.at("particles"))
        for (const auto& dist : sp.value("distribution", nlohmann::json::array()))
            if (dist.value("type", "") != "initial")
                throw std::invalid_argument("symplectic does not support injection");
    // Check the actual handler semantics, including omitted periodic axes.
    Domain checked_domain = domain_from_config(system_config);
    BoundaryConditionHandler checked_bc = boundary_conditions_from_config(system_config, checked_domain);
    for (int a = 0; a < 3; ++a)
        if (!checked_bc.is_periodic(a))
            throw std::invalid_argument("symplectic requires three periodic axes");
    if (system_config.value("StartFromTime", 0) > 0) {
        for (const auto& path : {"../Recovery/Fields/FieldE.backup", "../Recovery/Fields/FieldB.backup"})
            if (!std::filesystem::is_regular_file(path))
                throw std::runtime_error("missing symplectic recovery field");
        for (const auto& sp : particles_config.at("particles")) {
            const std::string name = sp.at("Name");
            if (!std::filesystem::is_regular_file("../Recovery/Particles/" + name + "/" + name + ".backup"))
                throw std::runtime_error("missing symplectic recovery species: " + name);
        }
    }
    Simulation::init();
    validate_state();
    for (auto& [name, sp] : species) sp->currentOnGrid.setZero();
    symplectic::WhitneyPIC pic(domain);
    Field3d rho(domain.size(), 1);
    collect_charge_density(rho);
    gauss_initial_.resize(domain.size(), 1);
    const int g = GHOST_CELLS;
    const auto n = domain.num_cells();
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k) {
                if (std::abs(pic.divergence(magnetic, i, j, k, true)) > 1e-12)
                    throw std::invalid_argument("symplectic requires divergence-free initial B");
                gauss_initial_(i + g, j + g, k + g, 0) =
                    pic.divergence(electric, i, j, k, false) - rho(i + g, j + g, k + g, 0);
            }
    energy_initial_ = pic.field_energy(electric) + pic.field_energy(magnetic);
    for (const auto& [name, sp] : species) energy_initial_ += sp->get_kinetic_energy();
}
void SimulationSymplectic::validate_state() const {
    RECORD_TIMER;
    if (!electric.data().allFinite() || !magnetic.data().allFinite())
        throw std::runtime_error("non-finite symplectic fields");
    for (const auto& [name, sp] : species)
        for (int cell = 0; cell < sp->size(); ++cell)
            for (const auto& p : (*sp)(cell)) {
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
        const auto cfg = system_config.value("Symplectic", nlohmann::json::object());
        const auto e = initial_vector(cfg, "initial_electric"), b = initial_vector(cfg, "initial_magnetic");
        for (int i = 0; i < electric.capacity(); ++i) {
            electric[i] = e[i % 3];
            magnetic[i] = b[i % 3];
        }
    }
    symplectic::WhitneyPIC(domain).sync(electric);
    symplectic::WhitneyPIC(domain).sync(magnetic);
}
int SimulationSymplectic::reset_deposit_buffers(std::vector<Field3d>& buffers, int components,
                                                const ParticlesArray& sp) const {
    RECORD_TIMER;
    const int threads = sp.get_total_num_of_particles() >= 4096 ? omp_get_max_threads() : 1;
    if (static_cast<int>(buffers.size()) < threads)
        buffers.resize(threads);
    for (int t = 0; t < threads; ++t) {
        if (buffers[t].sizes() != domain.size() || buffers[t].nd() != components)
            buffers[t].resize(domain.size(), components);
        buffers[t].setZero();
    }
    return threads;
}
void SimulationSymplectic::deposit_charge_density(const ParticlesArray& sp, Field3d& rho) const {
    RECORD_TIMER;
    symplectic::WhitneyPIC pic(domain);
    const double weight = sp.is_neutral() ? sp.mpw() : sp.charge() * sp.mpw();
    const int threads = reset_deposit_buffers(density_buffers_, 1, sp);
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        auto& local = density_buffers_[omp_get_thread_num()];
#pragma omp for schedule(static)
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const auto& p : sp(cell)) pic.deposit_charge(local, p.coord, weight);
    }
    if (rho.sizes() != domain.size() || rho.nd() != 1)
        rho.resize(domain.size(), 1);
    rho.setZero();
    for (int i = 0; i < rho.capacity(); ++i) {
        double value = 0;
        for (int t = 0; t < threads; ++t) value += density_buffers_[t][i];
        rho[i] = value;
    }
    pic.sync(rho);
}
void SimulationSymplectic::collect_charge_density(Field3d& field) const {
    RECORD_TIMER;
    field.resize(domain.size(), 1);
    field.setZero();
    for (const auto& [name, sp] : species) {
        if (sp->is_neutral())
            continue;
        Field3d rho(domain.size(), 1);
        deposit_charge_density(*sp, rho);
        field.data() += rho.data();
    }
}
void SimulationSymplectic::prepare_step(int timestep) {
    RECORD_TIMER;
    for (auto& [name, sp] : species) {
        sp->begin_step();
        sp->step_state.injection_energy =
            inject_particles_step(*sp, particle_sources.at(sp->name()).injection, timestep, domain, dt_,
                                  system_config.value("co_locate_species", true));
    }
}
void SimulationSymplectic::kinetic_flow(int axis, double dt) {
    RECORD_TIMER;
    symplectic::WhitneyPIC pic(domain);
    for (auto& [name, sp] : species) {
        const double qm = sp->is_neutral() ? 0 : sp->charge() / sp->mass();
        const double qw = sp->is_neutral() ? 0 : sp->charge() * sp->mpw();
        const int threads = reset_deposit_buffers(current_buffers_, 3, *sp);
        std::exception_ptr failure;
#pragma omp parallel num_threads(threads) if (threads > 1)
        {
            auto& integral = current_buffers_[omp_get_thread_num()];
#pragma omp for schedule(static)
            for (int cell = 0; cell < sp->size(); ++cell) {
                try {
                    for (auto& p : (*sp)(cell)) pic.advance_axis(p, axis, dt, qm, qw, magnetic, integral);
                } catch (...) {
#pragma omp critical(symplectic_trajectory_failure)
                    {
                        if (!failure)
                            failure = std::current_exception();
                    }
                }
            }
        }
        if (failure)
            std::rethrow_exception(failure);
        // H_axis holds B fixed; E does not enter this particle flow. Deposits
        // therefore commute and can be summed before updating the common E.
        for (int i = 0; i < electric.capacity(); ++i) {
            double integral = 0;
            for (int t = 0; t < threads; ++t) integral += current_buffers_[t][i];
            electric[i] -= integral;
            sp->currentOnGrid[i] += integral;
        }
    }
}
void SimulationSymplectic::electric_flow(double dt) {
    RECORD_TIMER;
    symplectic::WhitneyPIC pic(domain);
    pic.electric_flow(magnetic, electric, dt);
    for (auto& [name, sp] : species) {
        if (sp->is_neutral())
            continue;
        const double impulse = dt * sp->charge() / sp->mass();
        const int threads = sp->get_total_num_of_particles() >= 4096 ? omp_get_max_threads() : 1;
#pragma omp parallel for num_threads(threads) if (threads > 1) schedule(static)
        for (int cell = 0; cell < sp->size(); ++cell)
            for (auto& p : (*sp)(cell)) p.velocity += impulse * pic.gather(electric, p.coord, false);
    }
}
void SimulationSymplectic::make_step(int) {
    RECORD_TIMER;
    symplectic::WhitneyPIC pic(domain);
    // Xiao et al. (2015), Eq. (55): palindrome of exact Hamiltonian flows.
    // Unlike ECSIM, no half-time particle positions or averaged E are stored.
    kinetic_flow(0, 0.5 * dt_);
    kinetic_flow(1, 0.5 * dt_);
    kinetic_flow(2, 0.5 * dt_);
    pic.magnetic_flow(electric, magnetic, 0.5 * dt_);
    electric_flow(dt_);
    pic.magnetic_flow(electric, magnetic, 0.5 * dt_);
    kinetic_flow(2, 0.5 * dt_);
    kinetic_flow(1, 0.5 * dt_);
    kinetic_flow(0, 0.5 * dt_);
    pic.sync(electric);
    pic.sync(magnetic);
    for (auto& [name, sp] : species) {
        sp->currentOnGrid.data() /= dt_;
        pic.sync(sp->currentOnGrid);
        // Rebin only particles whose final cell changed. Delaying insertion
        // until the scan ends prevents processing a moved particle twice and
        // supports periodic seams and arbitrary multi-cell displacements.
        timer::commonTimer rebin_timer("symplectic rebin");
        std::vector<Particle> moved;
        for (int cell = 0; cell < sp->size(); ++cell) {
            auto& bucket = (*sp)(cell);
            size_t p = 0;
            while (p < bucket.size()) {
                const auto target = domain.storage_cell_index(bucket[p].coord);
                if (domain.sind(target[0], target[1], target[2]) == cell) {
                    ++p;
                } else {
                    moved.push_back(bucket[p]);
                    bucket[p] = bucket.back();
                    bucket.pop_back();
                }
            }
        }
        sp->add_particles(moved);
        rebin_timer.finish();
    }
    validate_state();
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
    symplectic::WhitneyPIC pic(domain);
    double kinetic = 0;
    for (const auto& [name, sp] : species) {
        const double energy = sp->get_kinetic_energy();
        diagnostic_->addEnergy(name, energy);
        kinetic += energy;
    }
    const double e = pic.field_energy(electric), b = pic.field_energy(magnetic);
    diagnostic_->addEnergy("energyFieldE", e);
    diagnostic_->addEnergy("energyFieldB", b);
    if (!std::isfinite(kinetic + e + b))
        throw std::runtime_error("non-finite symplectic energy");
    diagnostic_->addEnergy("energyTotal", kinetic + e + b);
    diagnostic_->addEnergy("energyError", kinetic + e + b - energy_initial_);
    Field3d rho(domain.size(), 1);
    collect_charge_density(rho);
    double gauss_drift = 0, gauss_error = 0, div_b = 0;
    const auto n = domain.num_cells();
    const int g = GHOST_CELLS;
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k) {
                const double defect = pic.divergence(electric, i, j, k, false) - rho(i + g, j + g, k + g, 0);
                gauss_error = std::max(gauss_error, std::abs(defect));
                gauss_drift = std::max(gauss_drift, std::abs(defect - gauss_initial_(i + g, j + g, k + g, 0)));
                div_b = std::max(div_b, std::abs(pic.divergence(magnetic, i, j, k, true)));
            }
    diagnostic_->addEnergy("gaussError", gauss_error);
    diagnostic_->addEnergy("gaussDrift", gauss_drift);
    diagnostic_->addEnergy("divB", div_b);
    write_diagnostics(timestep, *diagnostic_,
                      {domain, species, electric, magnetic, zero_field_, zero_field_, bc_handler, SHAPE}, outputs_);
}
void SimulationSymplectic::finalize() {
    RECORD_TIMER;
    for (auto& out : outputs_) out->finalize();
    Simulation::finalize();
}
