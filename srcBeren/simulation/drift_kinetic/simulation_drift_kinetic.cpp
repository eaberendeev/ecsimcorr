#include "simulation_drift_kinetic.h"

#include <omp.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "runtime/profiler.h"
#include "simulation/eccapfim/algorithms_eccapfim.h"
#include "utils/containers.h"
#include "utils/util.h"

namespace {
bool species_uses_drift_kinetic(const nlohmann::json& particles_config, const std::string& name) {
    if (!particles_config.contains("particles") || !particles_config["particles"].is_array())
        return false;
    for (const auto& entry : particles_config["particles"]) {
        if (!entry.is_object() || entry.value("Name", std::string()) != name)
            continue;
        const auto model = entry.value("model", std::string("kinetic"));
        if (model != "kinetic" && model != "drift_kinetic")
            throw std::invalid_argument("particle model of " + name + " must be kinetic or drift_kinetic");
        return model == "drift_kinetic";
    }
    return false;
}

Vector3R vector_from_json(const nlohmann::json& value, const Vector3R& fallback) {
    if (!value.is_array() || value.size() != 3)
        return fallback;
    return Vector3R(value[0].get<double>(), value[1].get<double>(), value[2].get<double>());
}
}   // namespace

void SimulationDriftKinetic::init() {
    SimulationEccapfim::init();
    coupling_ = std::make_unique<driftKinetic::YeeCoupling>(domain, &bc_handler);
    const auto eccapfim_config = system_config.value("Eccapfim", nlohmann::json::object());
    // B^{n+1} enters the guiding-centre orbit, so B is a nonlinear unknown.
    if (!eccapfim_config.contains("coupled_fields"))
        coupled_fields_ = true;
    if (!coupled_fields_)
        throw std::invalid_argument("drift_kinetic requires Eccapfim.coupled_fields=true");
    for (int d = 0; d < 3; ++d)
        if (!bc_handler.is_periodic(d))
            throw std::invalid_argument("drift_kinetic guiding centres support only periodic boundaries");
    // Guiding centres have their own recovery files next to the kinetic ones.
    for (const auto& output :
         system_config.value("diagnostics", nlohmann::json::object()).value("outputs", nlohmann::json::array()))
        if (output.value("type", std::string()) == "recovery" && output.contains("interval"))
            recovery_interval_ = output.at("interval").get<int>();

    const auto kinetic = system_config.value("DriftKinetic", nlohmann::json::object());
    evolve_fields_ = kinetic.value("evolve_fields", true);
    const auto push = system_config.value("DriftKineticPush", nlohmann::json::object());
    // Sequential Picard with safeguarded Aitken relaxation (xpic
    // DRIFT_KINETIC_AITKEN=1) halves the orbit iterations of plain Picard.
    const std::string solver = push.value("solver", std::string("aitken"));
    if (solver != "picard" && solver != "aitken")
        throw std::invalid_argument("DriftKineticPush.solver must be picard or aitken");
    use_aitken_ = solver == "aitken";
    push_.set_tolerances(push.value("atol", 1e-16), push.value("rtol", 1e-13), push.value("max_iterations", 100));
    adaptive_.configure(push.value("substepping", false), push.value("max_substep_depth", 4));
    coord_is_gc_ = push.value("coord_is_gc", false);
    fail_on_terminal_nonconvergence_ = push.value("fail_on_terminal_nonconvergence", false);
    fail_on_final_nonconvergence_ = push.value("fail_on_final_nonconvergence", false);

    for (auto* field : {&fieldJ_guiding, &fieldM, &magnetic_begin_, &magnetic_end_}) {
        field->resize(domain.size(), 3);
        field->setZero();
    }
    convert_initial_particles();
    // Restart: the kinetic markers of a guiding-centre species are not stored
    // (the species list holds them as guiding centres), so their containers are
    // empty here and are filled from the guiding-centre recovery files.
    if (system_config.value("StartFromTime", 0) > 0)
        for (auto& [name, ptr] : particles_) {
            if (ptr->get_total_num_of_particles() > 0)
                throw std::logic_error("drift_kinetic restart: " + name + " already has markers");
            driftKinetic::read_guiding_centers(*ptr, recovery_file("..", name));
            std::cout << "  Loaded guiding centres " << name << " from recovery\n";
        }
    previous_energy_ = guiding_center_energy(fieldB + fieldBInit);
    for (const auto& [name, ptr] : particles_) {
        Field3d& current = species_current_[name];
        current.resize(domain.size(), 3);
        current.setZero();
        moment_adapters_[name] = std::make_unique<driftKinetic::GuidingCenterMoments>(name, particles_, *coupling_,
                                                                                      fieldBn, fieldBInit, current);
        moment_sources_[name] = moment_adapters_[name].get();
    }
}

void SimulationDriftKinetic::init_fields() {
    SimulationEccapfim::init_fields();
    // Uniform prescribed fields for tests; ExternalFieldE/B take precedence.
    const auto kinetic = system_config.value("DriftKinetic", nlohmann::json::object());
    if (kinetic.contains("electric") && !system_config.contains("ExternalFieldE")) {
        const Vector3R electric = vector_from_json(kinetic["electric"], Vector3R(0, 0, 0));
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k)
                    for (int c = 0; c < 3; ++c) fieldE_external(i, j, k, c) = electric[c];
    }
    if (kinetic.contains("magnetic") && !system_config.contains("ExternalFieldB")) {
        const Vector3R magnetic = vector_from_json(kinetic["magnetic"], Vector3R(0, 0, 1));
        if (!(magnetic.length() > 0.0) || !std::isfinite(magnetic.length()))
            throw std::invalid_argument("drift_kinetic requires a finite nonzero DriftKinetic.magnetic field");
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k)
                    for (int c = 0; c < 3; ++c) fieldBInit(i, j, k, c) = magnetic[c];
    }
}

// Kinetic markers become guiding centres: by default the marker is the
// particle position and the centre is shifted by the Larmor vector
// v x b / (q/m |B|) (xpic `from_particle_position`); with `coord_is_gc` the
// marker already is the centre. mu = m v_perp^2 / (2 |B|) uses B at the marker.
void SimulationDriftKinetic::convert_initial_particles() {
    const std::string scheme = system_config.value("Scheme", std::string("drift_kinetic"));
    const bool all_species = scheme == "drift_kinetic";
    const Field3d magnetic = fieldB + fieldBInit;
    for (auto it = species.begin(); it != species.end();) {
        const std::string name = it->first;
        auto& sp = *it->second;
        if (!all_species && !species_uses_drift_kinetic(particles_config, name)) {
            ++it;
            continue;
        }
        if (sp.is_neutral())
            throw std::invalid_argument("drift_kinetic requires charged particle species: " + name);
        auto particles =
            std::make_unique<driftKinetic::DriftKineticParticles>(domain, name, sp.charge(), sp.mass(), sp.mpw());
        const double qm = sp.charge() / sp.mass();
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const Particle& particle : sp(cell)) {
                const Vector3R B = coupling_->magnetic_at(magnetic, particle.coord);
                const double lenB = B.length();
                if (!(lenB > 0.0) || !std::isfinite(lenB))
                    throw std::invalid_argument("drift_kinetic requires a nonzero magnetic field at every marker");
                const Vector3R b = B / lenB;
                const double v_parallel = particle.velocity.dot(b);
                const Vector3R v_perpendicular = particle.velocity - v_parallel * b;
                driftKinetic::GuidingCenter guiding;
                guiding.position = particle.coord;
                if (!coord_is_gc_)
                    guiding.position += particle.velocity.cross(b) / (qm * lenB);
                guiding.position = bc_handler.wrap_periodic(guiding.position, domain);
                guiding.v_parallel = v_parallel;
                guiding.mu = 0.5 * sp.mass() * v_perpendicular.squaredNorm() / lenB;
                guiding.charge = sp.charge();
                guiding.mass = sp.mass();
                guiding.weight = sp.mpw();
                particles->add(guiding);
            }
        particles_.emplace(name, std::move(particles));
        push_stats_[name] = {};
        charged_species.erase(name);
        particle_sources.erase(name);
        it = species.erase(it);
    }
}

double SimulationDriftKinetic::guiding_center_energy(const Field3d& magnetic) const {
    double energy = 0.0;
    for (const auto& [name, ptr] : particles_) {
        const auto& sp = *ptr;
#pragma omp parallel for reduction(+ : energy) schedule(dynamic, 64)
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const auto& p : sp(cell)) {
                const double lenB = coupling_->magnetic_at(magnetic, p.position).length();
                energy += p.weight * (0.5 * p.mass * p.v_parallel * p.v_parallel + p.mu * lenB);
            }
    }
    return energy;
}

Vector3R SimulationDriftKinetic::extra_momentum() const {
    const Field3d magnetic = fieldBn + fieldBInit;
    Vector3R momentum(0, 0, 0);
    for (const auto& [name, ptr] : particles_)
        for (int cell = 0; cell < ptr->size(); ++cell)
            for (const auto& p : (*ptr)(cell))
                momentum +=
                    (p.mass * p.weight * p.v_parallel) * coupling_->magnetic_at(magnetic, p.position).normalized();
    return momentum;
}

driftKinetic::DriftKineticParticles& SimulationDriftKinetic::guiding_centers(const std::string& name) {
    return *particles_.at(name);
}

const driftKinetic::PushStatistics& SimulationDriftKinetic::push_statistics(const std::string& name) const {
    return push_stats_.at(name);
}

void SimulationDriftKinetic::prepare_step(int timestep) {
    SimulationEccapfim::prepare_step(timestep);
    previous_energy_ = guiding_center_energy(fieldB + fieldBInit);
    warned_this_step_ = false;
}

double SimulationDriftKinetic::extra_kinetic_energy() const {
    return guiding_center_energy(fieldBn + fieldBInit);
}

double SimulationDriftKinetic::extra_previous_kinetic_energy() const {
    return previous_energy_;
}

double SimulationDriftKinetic::extra_external_work() const {
    return dt_ * dot_product_sum(fieldE_external, fieldJ_guiding, bc_handler.active_range(domain.grid));
}

void SimulationDriftKinetic::begin_nonlinear_step() {
    SimulationEccapfim::begin_nonlinear_step();
    initial_particles_.clear();
    trial_particles_.clear();
    for (const auto& [name, ptr] : particles_) {
        initial_particles_[name] = ptr->clone();
        trial_particles_[name] = ptr->clone();
    }
    magnetic_begin_ = fieldB + fieldBInit;
}

void SimulationDriftKinetic::evaluate_current(const Field3d& electric, const Field3d& magnetic) {
    RECORD_TIMER;
    SimulationEccapfim::evaluate_current(electric, magnetic);   // kinetic species -> fieldJ

    // `magnetic` is the full B^{n+1/2} iterate; B^{n+1} follows from the midpoint rule.
    magnetic_end_ = 2.0 * magnetic - magnetic_begin_;
    const driftKinetic::MagneticLevels levels{&magnetic_begin_, &magnetic, &magnetic_end_};
    coupling_->prepare(levels);
    // Thread-private source arrays, summed per species after its push: no atomic updates.
    const int threads = omp_get_max_threads();
    thread_current_.resize(threads);
    thread_magnetization_.resize(threads);
    for (int t = 0; t < threads; ++t)
        for (Field3d* field : {&thread_current_[t], &thread_magnetization_[t]})
            if (field->capacity() != fieldJ_guiding.capacity())
                field->resize(domain.size(), 3);
    fieldJ_guiding.setZero();
    fieldM.setZero();

    for (const auto& [name, initial_ptr] : initial_particles_) {
#pragma omp parallel for schedule(static)
        for (int t = 0; t < threads; ++t) {
            thread_current_[t].setZero();
            thread_magnetization_[t].setZero();
        }
        const auto& initial = *initial_ptr;
        auto& trial = *trial_particles_.at(name);
        driftKinetic::PushStatistics stats;
        long invalid = 0;
        std::exception_ptr failure;
        const int chunk = algorithmsECCAPFIM::particle_cell_chunk(initial.size());
#pragma omp parallel
        {
            driftKinetic::PushStatistics local;
            std::vector<driftKinetic::Leaf> leaves;
            Field3d& current = thread_current_[omp_get_thread_num()];
            Field3d& magnetization = thread_magnetization_[omp_get_thread_num()];
            // B at the orbit start does not change during the orbit iteration.
            Vector3R cached_start(0, 0, 0), cached_field(0, 0, 0);
            double cached_theta = -1.0;
            const auto fields = [&](const Vector3R& a, const Vector3R& b, double theta0, double theta1) {
                if (theta0 != cached_theta || !(a == cached_start)) {
                    cached_field = coupling_->start_magnetic(levels, a, {theta0, theta1});
                    cached_start = a;
                    cached_theta = theta0;
                }
                return coupling_->interpolate(electric, levels, a, b, {theta0, theta1}, &cached_field);
            };
#pragma omp for schedule(dynamic, chunk) reduction(+ : invalid)
            for (int cell = 0; cell < initial.size(); ++cell) {
                try {
                    const auto& start = initial(cell);
                    auto& end = trial(cell);
                    for (size_t p = 0; p < start.size(); ++p) {
                        if (!adaptive_.advance(push_, start[p], end[p], dt_, fields, use_aitken_, leaves, local)) {
                            ++invalid;
                            continue;
                        }
                        const double qw = start[p].charge * start[p].weight;
                        const double muw = start[p].mu * start[p].weight;
                        for (const auto& leaf : leaves) {
                            const double fraction = leaf.end - leaf.begin;
                            const Vector3R velocity = (leaf.finish.position - leaf.start.position) / (fraction * dt_);
                            coupling_->deposit(current, magnetization, levels, leaf.start.position,
                                               leaf.finish.position, velocity, qw * fraction, muw * fraction,
                                               {leaf.begin, leaf.end});
                        }
                        end[p] = leaves.back().finish;
                    }
                } catch (...) {
#pragma omp critical(drift_kinetic_failure)
                    if (!failure)
                        failure = std::current_exception();
                }
            }
#pragma omp critical(drift_kinetic_stats)
            stats.merge(local);
        }
        if (failure)
            std::rethrow_exception(failure);
        if (invalid > 0)
            throw std::runtime_error("drift_kinetic: invalid guiding-centre orbit for " + std::to_string(invalid) +
                                     " particles of " + name);
        if (stats.leaf_failures > 0) {
            if (fail_on_terminal_nonconvergence_)
                throw std::runtime_error("drift_kinetic: guiding-centre push did not converge for " + name);
            if (!warned_this_step_) {
                std::cerr << "WARNING: drift_kinetic push non-convergence step=" << current_timestep_
                          << " species=" << name << " leaf_failures=" << stats.leaf_failures
                          << " max_FR=" << stats.max_leaf_residue_r << " max_FV=" << stats.max_leaf_residue_v << '\n';
                warned_this_step_ = true;
            }
        }
        push_stats_[name] = stats;

        Field3d& species_current = species_free_current_[name];
        Field3d& species_magnetization = species_magnetization_[name];
        for (Field3d* field : {&species_current, &species_magnetization})
            if (field->capacity() != fieldJ_guiding.capacity())
                field->resize(domain.size(), 3);
#pragma omp parallel for schedule(static)
        for (int p = 0; p < fieldJ_guiding.capacity(); ++p) {
            double j = 0.0, m = 0.0;
            for (int t = 0; t < threads; ++t) {
                j += thread_current_[t](p);
                m += thread_magnetization_[t](p);
            }
            species_current(p) = j;
            species_magnetization(p) = m;
            fieldJ_guiding(p) += j;
            fieldM(p) += m;
        }
    }
    algorithmsECCAPFIM::copy_periodic_nodes(fieldJ_guiding, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(fieldM, domain, &bc_handler);
    // Magnetization current curl_B M on the same Yee operator as curl_B B.
    fieldJ += fieldJ_guiding + maxwell.curlB * fieldM;
}

void SimulationDriftKinetic::accept_nonlinear_step() {
    SimulationEccapfim::accept_nonlinear_step();
    for (auto& [name, trial] : trial_particles_) {
        const auto& stats = push_stats_.at(name);
        if (fail_on_final_nonconvergence_ && stats.failed_particles > 0)
            throw std::runtime_error("drift_kinetic: accepted push has " + std::to_string(stats.failed_particles) +
                                     " non-converged particles of " + name);
        for (int cell = 0; cell < trial->size(); ++cell)
            for (auto& p : (*trial)(cell)) p.position = bc_handler.wrap_periodic(p.position, domain);
        trial->update_cells();
    }
    particles_ = std::move(trial_particles_);
    trial_particles_.clear();
    initial_particles_.clear();
    // Diagnostic current of each species: J_gc + curl_B M of the accepted evaluation.
    for (auto& [name, current] : species_current_) {
        Field3d free = species_free_current_.at(name), magnetization = species_magnetization_.at(name);
        algorithmsECCAPFIM::copy_periodic_nodes(free, domain, &bc_handler);
        algorithmsECCAPFIM::copy_periodic_nodes(magnetization, domain, &bc_handler);
        current = free + maxwell.curlB * magnetization;
    }
}

// Prescribed fields (`DriftKinetic.evolve_fields = false`): one push with
// E^n and B^n, no field solve.
void SimulationDriftKinetic::make_step(int timestep) {
    if (evolve_fields_) {
        SimulationEccapfim::make_step(timestep);
        return;
    }
    RECORD_TIMER;
    current_timestep_ = timestep;
    residual_evaluations_ = 1;
    begin_nonlinear_step();
    fieldEp = fieldE;
    fieldBHalf = fieldB;
    evaluate_current(fieldE + fieldE_external, fieldB + fieldBInit);
    fieldEn = fieldE;
    fieldBn = fieldB;
    nonlinear_iterations = 0;
    nonlinear_residual = 0;
    accept_nonlinear_step();
}

// Guiding centres: V = v_par b + E x b / B with the mid-step v_par, so
// (dt/2) dJ/dE = (dt/2)^2 q^2 w / m b b^T + (dt/2) q w / B [. x b]
// (xpic hybrid_drift_kinetic::form_mass_matrix, straight field lines).
void SimulationDriftKinetic::assemble_mass_matrix(const Field3d& magnetic) {
    SimulationEccapfim::assemble_mass_matrix(magnetic);
    const bool have_trial = !trial_particles_.empty();
    const Species& source = have_trial ? trial_particles_ : particles_;
    for (const auto& [name, ptr] : source) {
        const auto& sp = *ptr;
        const auto& start = have_trial ? *initial_particles_.at(name) : sp;
        const double q = sp.charge(), m = sp.mass(), w = sp.weight();
#pragma omp parallel for schedule(dynamic, 16)
        for (int cell = 0; cell < sp.size(); ++cell)
            for (size_t p = 0; p < sp(cell).size(); ++p) {
                const Vector3R point = mass_matrix_point(start(cell)[p].position, sp(cell)[p].position);
                const Vector3R B = coupling_->magnetic_at(magnetic, point);
                const double lenB = B.length();
                if (!(lenB > 0.0) || !std::isfinite(lenB))
                    continue;
                const Vector3R b = B / lenB;
                const double a = 0.25 * dt_ * dt_ * q * q * w / m;
                const double c = 0.5 * dt_ * q * w / lenB;
                Eigen::Matrix3d T;
                T << a * b.x() * b.x(), a * b.x() * b.y() + c * b.z(), a * b.x() * b.z() - c * b.y(),
                    a * b.y() * b.x() - c * b.z(), a * b.y() * b.y(), a * b.y() * b.z() + c * b.x(),
                    a * b.z() * b.x() + c * b.y(), a * b.z() * b.y() - c * b.x(), a * b.z() * b.z();
                add_particle_mass(bc_handler.wrap_periodic(point, domain), T);
            }
    }
}

void SimulationDriftKinetic::add_extra_diagnostics(Diagnostics& diag) {
    const Field3d magnetic = fieldBn + fieldBInit;
    for (const auto& [name, ptr] : particles_) {
        const auto& sp = *ptr;
        double energy = 0.0;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const auto& p : sp(cell))
                energy += p.weight * (0.5 * p.mass * p.v_parallel * p.v_parallel +
                                      p.mu * coupling_->magnetic_at(magnetic, p.position).length());
        const auto& stats = push_stats_.at(name);
        diag.addEnergy(name, energy);
        diag.addEnergy(name + "Particles", static_cast<double>(sp.get_total_num_of_particles()));
        diag.addEnergy(name + "PushIterations", stats.average_iterations());
        diag.addEnergy(name + "PushMaxIterations", stats.max_iterations);
        diag.addEnergy(name + "PushRetries", static_cast<double>(stats.retries));
        diag.addEnergy(name + "PushLeafFailures", static_cast<double>(stats.leaf_failures));
        diag.addEnergy(name + "PushFailedParticles", static_cast<double>(stats.failed_particles));
    }
    // Guiding-centre side of the energy theorem: dK_gc = dt E . J_gc - M . dB.
    const auto range = bc_handler.active_range(domain.grid);
    const double current_work = dt_ * dot_product_sum(fieldEp + fieldE_external, fieldJ_guiding, range);
    const double magnetization_work = dot_product_sum(fieldM, fieldBn - fieldB, range);
    const double kinetic_change = guiding_center_energy(magnetic) - previous_energy_;
    diag.addEnergy("guidingCenterCurrentWork", current_work);
    diag.addEnergy("guidingCenterMagnetizationWork", magnetization_work);
    diag.addEnergy("guidingCenterEnergyResidual", kinetic_change - current_work + magnetization_work);
}

void SimulationDriftKinetic::make_diagnostic(int timestep) {
    SimulationEccapfim::make_diagnostic(timestep);
    // With the kinetic recovery output: the same step, the accepted r^{n+1}, v_par^{n+1}, mu.
    if (recovery_interval_ > 0 && timestep % recovery_interval_ == 0)
        for (const auto& [name, ptr] : particles_) {
            ensure_directory(".//Recovery//Particles//" + name);
            driftKinetic::write_guiding_centers(*ptr, recovery_file(".", name));
        }
}

std::string SimulationDriftKinetic::recovery_file(const std::string& root, const std::string& name) {
    return root + "//Recovery//Particles//" + name + "//" + name + ".gc.backup";
}
