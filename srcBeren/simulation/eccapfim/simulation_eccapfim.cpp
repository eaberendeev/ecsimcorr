#include "simulation_eccapfim.h"

#include <Eigen/QR>
#include <exception>
#include <iostream>
#include <stdexcept>

#include "Diagnostic.h"
#include "ParticlesDiagnostic.h"
#include "algorithms_eccapfim.h"
#include "external_fieldsB.h"
#include "external_fieldsE.h"
#include "recovery.h"

SimulationEccapfim::~SimulationEccapfim() = default;

void SimulationEccapfim::init() {
    const double dt = system_config.at("Dt").get<double>();
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("ECCAPFIM requires finite Dt > 0");
    for (const auto* key : {"Dx", "Dy", "Dz"}) {
        const double h = system_config.at(key).get<double>();
        if (!std::isfinite(h) || h <= 0)
            throw std::invalid_argument("ECCAPFIM requires positive grid spacing");
    }
    for (const auto* key : {"NumCellsX", "NumCellsY", "NumCellsZ"})
        if (system_config.at(key).get<int>() < 1)
            throw std::invalid_argument("ECCAPFIM requires positive grid sizes");
    if (system_config.contains("CylinderDomain") || system_config.value("DampingType", std::string("None")) != "None")
        throw std::invalid_argument("ECCAPFIM requires a rectangular box without damping");
    if (!system_config.contains("Boundary_conditions") || !system_config["Boundary_conditions"].is_array())
        throw std::invalid_argument("ECCAPFIM requires Boundary_conditions");
    for (const auto& entry : system_config["Boundary_conditions"]) {
        if (!entry.is_object() || entry.size() != 1 ||
            !(entry.contains("periodic") || entry.contains("open") || entry.contains("electron_reflection")))
            throw std::invalid_argument("ECCAPFIM supports periodic, open and electron_reflection boundaries");
        if (entry.contains("electron_reflection")) {
            const auto check_face = [](const auto& params) {
                const auto face = params.at("face").template get<std::string>();
                if (face != "ZMIN" && face != "ZMAX")
                    throw std::invalid_argument("ECCAPFIM electron_reflection requires a Z face");
            };
            const auto& params = entry.at("electron_reflection");
            if (params.is_array())
                for (const auto& face : params) check_face(face);
            else
                check_face(params);
        }
    }

    const auto cfg = system_config.value("Eccapfim", nlohmann::json::object());
    absolute_tolerance_ = cfg.value("atol", 1e-10);
    relative_tolerance_ = cfg.value("rtol", 1e-10);
    particle_tolerance_ = cfg.value("particle_tolerance", 1e-12);
    mixing_ = cfg.value("mixing", 0.5);
    max_iterations_ = cfg.value("max_iterations", 200);
    history_ = cfg.value("history", 6);
    particle_iterations_ = cfg.value("particle_iterations", 100);
    if (!std::isfinite(absolute_tolerance_) || absolute_tolerance_ <= 0 || !std::isfinite(relative_tolerance_) ||
        relative_tolerance_ < 0 || !std::isfinite(particle_tolerance_) || particle_tolerance_ <= 0 ||
        !std::isfinite(mixing_) || mixing_ <= 0 || mixing_ > 1 || max_iterations_ < 1 || history_ < 0 ||
        history_ > 30 || particle_iterations_ < 1)
        throw std::invalid_argument("ECCAPFIM: invalid nonlinear solver settings");
    Simulation::init();
}

SimulationEccapfim::Vector SimulationEccapfim::pack(const Field3d& field) const {
    const auto n = unknown_range_.end - unknown_range_.start;
    const auto first = unknown_range_.start;
    Vector result(3 * n.elements_product());
    int p = 0;
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k)
                for (int c = 0; c < 3; ++c) result[p++] = field(i + first.x(), j + first.y(), k + first.z(), c);
    return result;
}

void SimulationEccapfim::unpack(const Vector& vector, Field3d& field) const {
    const auto n = unknown_range_.end - unknown_range_.start;
    const auto first = unknown_range_.start;
    int p = 0;
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k)
                for (int c = 0; c < 3; ++c) field(i + first.x(), j + first.y(), k + first.z(), c) = vector[p++];
    algorithmsECCAPFIM::copy_periodic_nodes(field, domain, &bc_handler);
}

void SimulationEccapfim::init_operators() {
    // Every nonperiodic face must define the fate of an outgoing particle.
    for (int d = 0; d < 3; ++d) {
        const Face low[] = {Face::XMIN, Face::YMIN, Face::ZMIN};
        const Face high[] = {Face::XMAX, Face::YMAX, Face::ZMAX};
        if (!bc_handler.is_periodic(d) &&
            (!bc_handler.path_policy().clips(low[d]) || !bc_handler.path_policy().clips(high[d])))
            throw std::invalid_argument("ECCAPFIM requires open or reflecting conditions on both nonperiodic faces");
        if (bc_handler.is_periodic(d) &&
            (bc_handler.path_policy().clips(low[d]) || bc_handler.path_policy().clips(high[d])))
            throw std::invalid_argument("ECCAPFIM cannot combine periodic and removing conditions on one axis");
    }
    unknown_range_ = bc_handler.active_range(domain.grid);
    const auto n = unknown_range_.end - unknown_range_.start;
    const int g = domain.grid.ghost_cells();
    auto reduced_index = [&](int full) {
        int index[3];
        for (int d = 0; d < 3; ++d)
            index[d] = algorithmsECCAPFIM::node_owner(domain.pos_vind(full, d) - g, d, domain, &bc_handler) -
                       unknown_range_.start[d];
        return 3 * (index[0] * n.y() * n.z() + index[1] * n.z() + index[2]) + domain.pos_vind(full, 3);
    };
    Field3d mask(domain.size(), 3);
    mask.data().setOnes();
    // Reuse the exact component-dependent Yee mask of open/electron_reflection.
    // ELECTRIC does not sum periodic copies; deposits already target owners.
    bc_handler.apply_to_fields(mask, FieldType::ELECTRIC, domain);
    electric_free_ = pack(mask);
    density_free_.resize(domain.size(), 1);
    density_free_.data().setOnes();
    bc_handler.apply_to_fields(density_free_, FieldType::DENSITY, domain);
    for (int p = 0; p < density_free_.capacity(); ++p)
        density_free_[p] = density_free_[p] != 0 ? 1 : 0;

    std::vector<Eigen::Triplet<double>> entries;
    // Restrict the existing Yee M = I + dt^2/4 curlB curlE to independent
    // nodes. Zero Dirichlet DOFs have identity rows and no free-row columns.
    // Thus operator, RHS and nonlinear current impose the same E_b = 0.
    for (int outer = 0; outer < mesh.IMmat.outerSize(); ++outer)
        for (Operator::InnerIterator it(mesh.IMmat, outer); it; ++it) {
            bool physical = true;
            for (int d = 0; d < 3; ++d) {
                const int i = domain.pos_vind(it.row(), d);
                physical = physical && i >= unknown_range_.start[d] && i < unknown_range_.end[d];
            }
            if (!physical)
                continue;
            const int row = reduced_index(it.row()), col = reduced_index(it.col());
            if (electric_free_[row] != 0 && electric_free_[col] != 0)
                entries.emplace_back(row, col, it.value());
        }
    for (int row = 0; row < electric_free_.size(); ++row)
        if (electric_free_[row] == 0)
            entries.emplace_back(row, row, 1.0);
    maxwell_.resize(3 * n.elements_product(), 3 * n.elements_product());
    maxwell_.setFromTriplets(entries.begin(), entries.end());
    maxwell_solver_.setTolerance(1e-13);
    maxwell_solver_.setMaxIterations(2000);
    maxwell_solver_.compute(maxwell_);
    if (maxwell_solver_.info() != Eigen::Success)
        throw std::runtime_error("ECCAPFIM: invalid Maxwell operator");
}

void SimulationEccapfim::init_fields() {
    for (auto* field : {&fieldE, &fieldEn, &fieldEp, &fieldB, &fieldBn, &fieldBInit, &fieldE_external, &fieldJ}) {
        field->resize(domain.size(), 3);
        field->setZero();
    }
    if (system_config.value("StartFromTime", 0) > 0)
        read_fields_from_recovery(fieldE, fieldB);
    if (auto cfg = create_electric_field_config(system_config, "ExternalFieldE"))
        cfg->apply(fieldE_external, domain);
    if (auto cfg = create_magnetic_field_config(system_config, "ExternalFieldB"))
        cfg->apply(fieldBInit, domain);
    for (auto* field : {&fieldE, &fieldB, &fieldE_external, &fieldBInit}) {
        if (!field->data().allFinite())
            throw std::invalid_argument("ECCAPFIM: nonfinite initial fields");
        algorithmsECCAPFIM::copy_periodic_nodes(*field, domain, &bc_handler);
    }
    bc_handler.apply_to_fields(fieldE, FieldType::ELECTRIC, domain);
    bc_handler.apply_to_fields(fieldB, FieldType::MAGNETIC, domain);
    fieldEn = fieldE;
    fieldBn = fieldB;
    fieldEp = fieldE;
}

void SimulationEccapfim::deposit_density(ParticlesArray& sp) {
    sp.densityOnGrid.setZero();
    // Neutral density is number density. It is excluded from charge collection.
    const double weight = sp.mpw() * (sp.is_neutral() ? 1 : sp.charge);
#pragma omp parallel for
    for (int cell = 0; cell < sp.size(); ++cell)
        for (const auto& particle : sp.particlesData(cell))
            algorithmsECCAPFIM::deposit_charge(sp.densityOnGrid, particle.coord, weight, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(sp.densityOnGrid, domain, &bc_handler);
    sp.densityOnGrid.data().array() *= density_free_.data().array();
}

void SimulationEccapfim::prepare_step(int timestep) {
    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJ.setZero();
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        sp.prepare();
        sp.diag.injection_energy = sp.inject_particles_step(sp.get_injection_distributions(), timestep, domain, dt_,
                                                            system_config.value("co_locate_species", true));
    }
}

SimulationEccapfim::ParticleIterations SimulationEccapfim::save_particles() {
    ParticleIterations states;
    for (auto& [name, ptr] : species) {
        auto& state = states[name];
        state.species = ptr.get();
        state.current.resize(domain.size(), 3);
        state.initial.resize(domain.size());
        state.trial.resize(domain.size());
        state.events.resize(domain.size());
        state.initial = ptr->particlesData;
        state.trial = ptr->particlesData;
        for (int cell = 0; cell < state.initial.capacity(); ++cell)
            state.events(cell).resize(state.initial(cell).size());
    }
    return states;
}

void SimulationEccapfim::push_particles_and_deposit_current(ParticleIterations& states, const Field3d& electric,
                                                            const Field3d& magnetic) {
    fieldJ.setZero();
    for (auto& [name, state] : states) {
        auto& sp = *state.species;
        state.current.setZero();
        const auto policy = bc_handler.path_policy(sp, domain);
        std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 64)
        for (int cell = 0; cell < state.initial.capacity(); ++cell) {
            try {
                const auto& initial = state.initial(cell);
                auto& trial = state.trial(cell);
                for (size_t p = 0; p < initial.size(); ++p)
                    trial[p] = algorithmsECCAPFIM::advance_particle(initial[p], electric, magnetic, state.current,
                                                                    sp.is_neutral() ? 0 : sp.charge / sp.mass(),
                                                                    sp.is_neutral() ? 0 : sp.charge * sp.mpw(), dt_,
                                                                    domain, particle_tolerance_, particle_iterations_,
                                                                    &bc_handler, policy, &state.events(cell)[p]);
            } catch (...) {
#pragma omp critical(eccapfim_particle_failure)
                {
                    if (!failure)
                        failure = std::current_exception();
                }
            }
        }
        if (failure)
            std::rethrow_exception(failure);
        algorithmsECCAPFIM::copy_periodic_nodes(state.current, domain, &bc_handler);
        fieldJ += state.current;
    }
}

void SimulationEccapfim::solve_electric_field(ParticleIterations& states) {
    const Field3d magnetic = fieldB + fieldBInit;   // Lorentz force uses B^n + B_external.
    const Vector rhs = pack(fieldE + (0.5 * dt_) * (mesh.curlB * fieldB)).cwiseProduct(electric_free_);
    auto residual = [&](const Vector& electric) -> Vector {
        if (!electric.allFinite())
            throw std::runtime_error("ECCAPFIM: nonfinite field iterate");
        unpack(electric, fieldEp);
        const Field3d full_electric = fieldEp + fieldE_external;
        push_particles_and_deposit_current(states, full_electric, magnetic);
        // Maxwell residual F(Eh) = M Eh - E^n - dt/2 curl B^n + dt/2 J[Eh].
        return maxwell_ * electric - rhs + (0.5 * dt_) * pack(fieldJ).cwiseProduct(electric_free_);
    };

    Vector x = pack(fieldE);
    std::vector<Vector> old_f, old_g;
    double initial_norm = 0;
    bool converged = false;
    for (int it = 0; it <= max_iterations_; ++it) {
        const Vector r = residual(x);
        nonlinear_residual = r.norm();
        nonlinear_iterations = it;
        if (!r.allFinite() || !std::isfinite(nonlinear_residual))
            throw std::runtime_error("ECCAPFIM: nonfinite nonlinear residual");
        if (it == 0)
            initial_norm = nonlinear_residual;
        if (nonlinear_residual <= absolute_tolerance_ + relative_tolerance_ * initial_norm) {
            converged = true;
            break;
        }
        if (it == max_iterations_)
            break;
        const Vector f = mixing_ * maxwell_solver_.solve(r);
        if (maxwell_solver_.info() != Eigen::Success || !f.allFinite())
            throw std::runtime_error("ECCAPFIM: Maxwell preconditioner did not converge");
        const Vector g = x - f;
        Vector next = g;
        // Anderson acceleration of the Maxwell-preconditioned fixed point.
        // QR handles dependent history columns without forming normal equations.
        const int m = static_cast<int>(old_f.size());
        if (m > 0) {
            Eigen::MatrixXd df(x.size(), m), dg(x.size(), m);
            for (int j = 0; j < m; ++j) {
                df.col(j) = old_f[j] - f;
                dg.col(j) = old_g[j] - g;
            }
            const Vector coefficients = df.colPivHouseholderQr().solve(-f);
            const Vector accelerated = g + dg * coefficients;
            if (accelerated.allFinite() && (accelerated - x).norm() <= 10 * f.norm())
                next = accelerated;
            else {
                old_f.clear();
                old_g.clear();
            }
        }
        if (history_ > 0) {
            old_f.push_back(f);
            old_g.push_back(g);
            if (static_cast<int>(old_f.size()) > history_) {
                old_f.erase(old_f.begin());
                old_g.erase(old_g.begin());
            }
        }
        x = next;
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: nonlinear field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
}

void SimulationEccapfim::accept_particles(ParticleIterations& states) {
    for (auto& [name, state] : states) {
        auto& sp = *state.species;
        // Full rebin: implicit orbits may cross several cells or the periodic seam.
        for (int cell = 0; cell < sp.size(); ++cell) sp.particlesData(cell).clear();
        for (int cell = 0; cell < state.trial.capacity(); ++cell)
            for (size_t p = 0; p < state.trial(cell).size(); ++p) {
                const auto& particle = state.trial(cell)[p];
                const auto& event = state.events(cell)[p];
                if (event.removed != Face::UNAVAILABLE)
                    sp.diag.add_loss(event.removed, get_energy_particle(particle.velocity, sp.mass(), sp.mpw()));
                else
                    sp.add_particle(particle);
                for (size_t f = 0; f < event.reflected.size(); ++f)
                    if (event.reflected[f] != 0)
                        sp.diag.boundary[static_cast<Face>(f)].reflected_count += event.reflected[f];
            }
        sp.currentOnGrid = state.current;
        deposit_density(sp);
    }
}

void SimulationEccapfim::make_step([[maybe_unused]] int timestep) {
    RECORD_TIMER;
    globalTimer.start("Total");
    // x^n and v^n remain unchanged while Eh, particle orbits and J are iterated.
    auto particles = save_particles();
    solve_electric_field(particles);

    // Eh = (E^{n+1} + E^n)/2; Faraday advances the self-consistent B field.
    fieldEn = 2.0 * fieldEp - fieldE;
    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, dt_);
    algorithmsECCAPFIM::copy_periodic_nodes(fieldEn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    algorithmsECCAPFIM::copy_periodic_nodes(fieldBn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);

    // The last residual used the accepted Eh: commit those orbits and their J once.
    accept_particles(particles);
    globalTimer.finish("Total");
    std::cout << "ECCAPFIM iterations=" << nonlinear_iterations << " residual=" << nonlinear_residual << '\n';
}

void SimulationEccapfim::make_diagnostic(int timestep) {
    if (!diagnostic_) {
        diagnostic_ = std::make_unique<Diagnostics>(species);
        outputs_ = OutputFactory::create(system_config.value("diagnostics", nlohmann::json::object()), domain, species,
                                         system_config, fieldEn, fieldBn, fieldE_external, fieldBInit);
    }
    auto& diag = *diagnostic_;
    const auto range = bc_handler.active_range(domain.grid);
    double kinetic = 0, previous = 0, injected = 0, external_work = 0, lost = 0;
    for (const auto& [name, ptr] : species) {
        const auto moments = ptr->diagnostic_moments();
        diag.addEnergy(name, moments.energy());
        diag.addEnergy(name + "Particles", static_cast<double>(moments.count));
        diag.addEnergy(name + "Inject", ptr->diag.injection_energy);
        diag.addEnergy(name + "XY", moments.kinetic.x() + moments.kinetic.y());
        diag.addEnergy(name + "Z", moments.kinetic.z());
        diag.addBoundary(name + "_InjectE", ptr->diag.injection_energy);
        diag.addBoundary(name + "_InjectN", ptr->diag.injection_count);
        const char* face_names[] = {"XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"};
        int fi = 0;
        for (const auto face : {Face::XMIN, Face::XMAX, Face::YMIN, Face::YMAX, Face::ZMIN, Face::ZMAX}) {
            const auto it = ptr->diag.boundary.find(face);
            const PerFaceStats stats = it != ptr->diag.boundary.end() ? it->second : PerFaceStats{};
            const std::string prefix = name + "_" + face_names[fi++];
            diag.addBoundary(prefix + "LostE", stats.lost_energy);
            diag.addBoundary(prefix + "LostN", stats.lost_count);
            diag.addBoundary(prefix + "ReflN", stats.reflected_count);
            diag.addBoundary(prefix + "EmitE", 0);
            diag.addBoundary(prefix + "EmitN", 0);
        }
        lost += ptr->diag.total_lost_energy();
        kinetic += moments.energy();
        previous += ptr->kineticEnergy;
        injected += ptr->diag.injection_energy;
        external_work += dt_ * dot_product_sum(fieldE_external, ptr->currentOnGrid, range);
    }
    const double e = calc_energy_field(fieldEn, range), b = calc_energy_field(fieldBn, range);
    const double delta = e + b - calc_energy_field(fieldE, range) - calc_energy_field(fieldB, range);
    const double balance = kinetic - previous - injected + delta - external_work + lost;
    diag.addEnergy("energyFieldE", e);
    diag.addEnergy("energyFieldB", b);
    diag.addEnergy("energyFieldBFull", calc_energy_field(fieldBn + fieldBInit, range));
    diag.addEnergy("totalInjectEnergy", injected);
    diag.addEnergy("totalEmitEnergy", 0);
    diag.addEnergy("totalLostEnergy", lost);
    diag.addEnergy("energyDamping", 0);
    diag.addEnergy("deltaKinetic", kinetic - previous);
    diag.addEnergy("deltaField", delta);
    diag.addEnergy("JeExWork", external_work);
    diag.addEnergy("JeRealWork", dt_ * dot_product_sum(fieldEp, fieldJ, range));
    diag.addEnergy("energyResidualRaw", balance);
    diag.addEnergy("energyResidualDiscrete", balance);
    diag.addEnergy("energyConserve", std::abs(balance));
    diag.addEnergy("nonlinearResidual", nonlinear_residual);
    diag.addEnergy("nonlinearIterations", nonlinear_iterations);
    for (auto& output : outputs_) output->output(timestep, diag);
}

void SimulationEccapfim::finalize() {
    for (auto& output : outputs_) output->finalize();
    Simulation::finalize();
}
