#include "simulation_eccapfim.h"

#include <omp.h>

#include <Eigen/Core>
#include <algorithm>
#include <Eigen/QR>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "diagnostics/Diagnostic.h"
#include "diagnostics/DiagnosticMeasurements.h"
#include "diagnostics/OutputFactory.h"
#include "config/field_config.h"
#include "particles/particle_sources.h"
#include "algorithms_eccapfim.h"
#include "coupling/stencil.h"
#include "fields/external_fieldsB.h"
#include "fields/external_fieldsE.h"
#include "io/recovery.h"
#include "io/field_io.h"
#include "runtime/profiler.h"

SimulationEccapfim::~SimulationEccapfim() = default;

void SimulationEccapfim::init() {
    RECORD_TIMER;
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
    failure_dump_ = cfg.value("failure_dump", false);
    iteration_log_ = cfg.value("iteration_log", false);
    const std::string preconditioner = cfg.value("preconditioner", "diagonal");
    if (preconditioner != "diagonal" && preconditioner != "incomplete_cholesky" && preconditioner != "amgcl" &&
        preconditioner != "block_jacobi" && preconditioner != "mass_matrix")
        throw std::invalid_argument(
            "Eccapfim.preconditioner must be diagonal, incomplete_cholesky, amgcl, block_jacobi or mass_matrix");
    use_ic_ = preconditioner == "incomplete_cholesky";
    use_amgcl_ = preconditioner == "amgcl";
    use_block_jacobi_ = preconditioner == "block_jacobi";
    use_mass_matrix_ = preconditioner == "mass_matrix";
    const std::string nonlinear_solver = cfg.value("nonlinear_solver", std::string("anderson"));
    if (nonlinear_solver != "anderson" && nonlinear_solver != "ngmres")
        throw std::invalid_argument("Eccapfim.nonlinear_solver must be anderson or ngmres");
    use_ngmres_ = nonlinear_solver == "ngmres";
    coupled_fields_ = cfg.value("coupled_fields", false);
    // xpic uses the mass matrix as the right nonlinear preconditioner of
    // NGMRES; Newton-Krylov with finite-difference Jacobian products remains
    // available and is the default without a mass matrix.
    const std::string coupled_solver =
        cfg.value("coupled_solver", std::string(use_mass_matrix_ ? "anderson" : "newton_krylov"));
    if (coupled_solver != "anderson" && coupled_solver != "ngmres" && coupled_solver != "newton_krylov")
        throw std::invalid_argument("Eccapfim.coupled_solver must be anderson, ngmres or newton_krylov");
    coupled_solver_ = coupled_solver == "anderson" ? CoupledSolver::Anderson
                      : coupled_solver == "ngmres" ? CoupledSolver::Ngmres
                                                   : CoupledSolver::NewtonKrylov;
    if (coupled_solver_ != CoupledSolver::NewtonKrylov && !use_mass_matrix_)
        throw std::invalid_argument("Eccapfim.coupled_solver=" + coupled_solver +
                                    " requires preconditioner=mass_matrix");
    coupled_history_ = cfg.value("coupled_history", 30);
    mass_matrix_every_iteration_ = cfg.value("mass_matrix_every_iteration", false);
    mass_matrix_check_ = cfg.value("mass_matrix_check", false);
    const std::string mass_point = cfg.value("mass_matrix_point", std::string("orbit_midpoint"));
    if (mass_point != "orbit_midpoint" && mass_point != "orbit_end")
        throw std::invalid_argument("Eccapfim.mass_matrix_point must be orbit_midpoint or orbit_end");
    mass_matrix_at_end_ = mass_point == "orbit_end";
    const std::string mass_solver = cfg.value("mass_matrix_solver", std::string("gmres"));
    if (mass_solver != "sparse_lu" && mass_solver != "gmres")
        throw std::invalid_argument("Eccapfim.mass_matrix_solver must be gmres or sparse_lu");
    mass_matrix_lu_ = mass_solver == "sparse_lu";
    mass_matrix_tolerance_ = cfg.value("mass_matrix_tolerance", 1e-6);
    mass_preconditioner_name_ = cfg.value("mass_matrix_preconditioner", std::string("auto"));
    if (mass_preconditioner_name_ != "auto" && mass_preconditioner_name_ != "block_jacobi" &&
        mass_preconditioner_name_ != "fft" && mass_preconditioner_name_ != "fft_block_jacobi")
        throw std::invalid_argument(
            "Eccapfim.mass_matrix_preconditioner must be auto, block_jacobi, fft or fft_block_jacobi");
    if (coupled_solver_ == CoupledSolver::Ngmres && coupled_history_ < 1)
        throw std::invalid_argument("Eccapfim.coupled_history must be at least 1 for NGMRES");
    if (coupled_history_ < 0 || coupled_history_ > 100 || !(mass_matrix_tolerance_ > 0) ||
        !std::isfinite(mass_matrix_tolerance_))
        throw std::invalid_argument("ECCAPFIM: invalid NGMRES or mass-matrix settings");
    coupled_atol_ = cfg.value("coupled_atol", 1e-10);
    coupled_rtol_ = cfg.value("coupled_rtol", 1e-10);
    coupled_max_outer_ = cfg.value("coupled_max_outer", 20);
    coupled_gmres_restart_ = cfg.value("coupled_gmres_restart", 20);
    coupled_gmres_max_ = cfg.value("coupled_gmres_max", 60);
    if (!std::isfinite(coupled_atol_) || coupled_atol_ <= 0 || !std::isfinite(coupled_rtol_) ||
        coupled_rtol_ < 0 || coupled_max_outer_ < 1 || coupled_gmres_restart_ < 1 || coupled_gmres_max_ < 1)
        throw std::invalid_argument("ECCAPFIM: invalid coupled-field solver settings");
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
    init_operators();
}

SimulationEccapfim::Vector SimulationEccapfim::pack(const Field3d& field) const {
    RECORD_TIMER;
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
    RECORD_TIMER;
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
    RECORD_TIMER;
    MaxwellStencils stencils = assemble_maxwell_stencils(domain, bc_handler);
    field_system_.init(stencils, domain, dt_, bc_handler);
    maxwell.init(std::move(stencils));
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
    const int g = GHOST_CELLS;
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
    for (int p = 0; p < density_free_.capacity(); ++p) density_free_[p] = density_free_[p] != 0 ? 1 : 0;

    std::vector<Eigen::Triplet<double>> entries;
    // Restrict the existing Yee M = I + dt^2/4 curlB curlE to independent
    // nodes. Zero Dirichlet DOFs have identity rows and no free-row columns.
    // Thus operator, RHS and nonlinear current impose the same E_b = 0.
    for (int outer = 0; outer < field_system_.IM.outerSize(); ++outer)
        for (Operator::InnerIterator it(field_system_.IM, outer); it; ++it) {
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
    if (use_mass_matrix_) {
        const Eigen::Vector3i first(unknown_range_.start.x(), unknown_range_.start.y(), unknown_range_.start.z());
        const Eigen::Vector3i nodes(n.x(), n.y(), n.z());
        bool periodic = electric_free_.minCoeff() != 0;
        for (int d = 0; d < 3; ++d) periodic = periodic && bc_handler.is_periodic(d);
        stencil_operator_.init(
            first, nodes, [&](int i, int j, int k, int c) { return reduced_dof(i, j, k, c); }, electric_free_, maxwell_,
            periodic);
        using P = YeeStencilMatrix::Preconditioner;
        if (mass_preconditioner_name_ == "auto")
            mass_preconditioner_kind_ = periodic ? P::Circulant : P::BlockJacobi;
        else
            mass_preconditioner_kind_ = mass_preconditioner_name_ == "fft"     ? P::Circulant
                                        : mass_preconditioner_name_ == "fft_block_jacobi" ? P::CirculantBlockJacobi
                                                                                 : P::BlockJacobi;
        if (mass_preconditioner_kind_ != P::BlockJacobi && !periodic)
            throw std::invalid_argument("Eccapfim.mass_matrix_preconditioner=fft needs a periodic unconstrained box");
    }
    maxwell_solver_.setTolerance(1e-13);
    maxwell_solver_.setMaxIterations(2000);
    maxwell_solver_.compute(maxwell_);
    if (maxwell_solver_.info() != Eigen::Success)
        throw std::runtime_error("ECCAPFIM: invalid Maxwell operator");
    if (use_ic_) {
        maxwell_ic_solver_.setTolerance(1e-13);
        maxwell_ic_solver_.setMaxIterations(2000);
        maxwell_ic_solver_.compute(maxwell_);
        if (maxwell_ic_solver_.info() != Eigen::Success)
            throw std::runtime_error("ECCAPFIM: incomplete Cholesky preconditioner failed");
    }
    if (use_amgcl_) {
        try {
            // The AMGCL Eigen adapter consumes rows; store a RowMajor copy.
            maxwell_row_ = maxwell_;
            AmgSolver::params params;
            params.solver.tol = 1e-13;
            params.solver.maxiter = 2000;
            amgcl_solver_ = std::make_unique<AmgSolver>(maxwell_row_, params);
        } catch (const std::exception& error) {
            throw std::runtime_error(std::string("ECCAPFIM: AMGCL setup failed: ") + error.what());
        }
    }
    if (use_block_jacobi_) {
        maxwell_block_solver_.setTolerance(1e-13);
        maxwell_block_solver_.setMaxIterations(2000);
        maxwell_block_solver_.compute(maxwell_);
        if (maxwell_block_solver_.info() != Eigen::Success)
            throw std::runtime_error("ECCAPFIM: block-Jacobi preconditioner failed");
    }
}

void SimulationEccapfim::init_fields() {
    RECORD_TIMER;
    for (auto* field :
         {&fieldE, &fieldEn, &fieldEp, &fieldB, &fieldBn, &fieldBHalf, &fieldBInit, &fieldE_external, &fieldJ}) {
        field->resize(domain.size(), 3);
        field->setZero();
    }
    if (system_config.value("StartFromTime", 0) > 0)
        read_fields_from_recovery(fieldE, fieldB);
    else {
        // Self-consistent E^0 and B^0, e.g. a wave perturbation.
        if (auto cfg = create_electric_field_config(system_config, "InitialFieldE"))
            cfg->apply(fieldE, domain);
        if (auto cfg = create_magnetic_field_config(system_config, "InitialFieldB"))
            cfg->apply(fieldB, domain);
    }
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

void SimulationEccapfim::collect_charge_density(Field3d& field) const {
    RECORD_TIMER;
    field.resize(domain.size(), 1);
    field.setZero();
    // Neutral species carry number density and do not enter the charge balance.
    for (const auto& [name, ptr] : species) {
        if (ptr->is_neutral())
            continue;
        const double weight = ptr->mpw() * ptr->charge();
        for (int cell = 0; cell < ptr->size(); ++cell)
            for (const auto& particle : (*ptr)(cell))
                algorithmsECCAPFIM::deposit_charge(field, particle.coord, weight, domain, &bc_handler);
    }
    algorithmsECCAPFIM::copy_periodic_nodes(field, domain, &bc_handler);
    field.data().array() *= density_free_.data().array();
}

void SimulationEccapfim::prepare_step(int timestep) {
    RECORD_TIMER;
    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJ.setZero();
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        sp.begin_step();
        sp.step_state.injection_energy =
            inject_particles_step(sp, particle_sources.at(sp.name()).injection, timestep, domain, dt_,
                                  system_config.value("co_locate_species", true));
    }
}

SimulationEccapfim::ParticleIterations SimulationEccapfim::save_particles() {
    RECORD_TIMER;
    ParticleIterations states;
    for (auto& [name, ptr] : species) {
        auto& state = states[name];
        state.species = ptr.get();
        state.current.resize(domain.size(), 3);
        state.initial.resize(domain.size());
        state.trial.resize(domain.size());
        state.events.resize(domain.size());
        state.initial = ptr->cells();
        state.trial = ptr->cells();
        for (int cell = 0; cell < state.initial.capacity(); ++cell)
            state.events(cell).resize(state.initial(cell).size());
    }
    return states;
}

void SimulationEccapfim::push_particles_and_deposit_current(ParticleIterations& states, const Field3d& electric,
                                                            const Field3d& magnetic) {
    RECORD_TIMER;
    fieldJ.setZero();
    // Thread-private current arrays, summed after the push: no atomic updates.
    const int threads = omp_get_max_threads();
    kinetic_thread_current_.resize(threads);
    for (auto& current : kinetic_thread_current_)
        if (current.capacity() != fieldJ.capacity())
            current.resize(domain.size(), 3);
    for (auto& [name, state] : states) {
        auto& sp = *state.species;
#pragma omp parallel for schedule(static)
        for (int t = 0; t < threads; ++t) kinetic_thread_current_[t].setZero();
        const auto policy = bc_handler.path_policy(sp, domain);
        std::exception_ptr failure;
        int failed_cell = -1;
        size_t failed_particle = 0;
        {
            RECORD_TIMER_NAMED("eccapfim particle push");
            const int chunk = algorithmsECCAPFIM::particle_cell_chunk(state.initial.capacity());
#pragma omp parallel for schedule(dynamic, chunk)
            for (int cell = 0; cell < state.initial.capacity(); ++cell) {
                size_t particle_index = 0;
                try {
                    const auto& initial = state.initial(cell);
                    auto& trial = state.trial(cell);
                    for (size_t p = 0; p < initial.size(); ++p) {
                        particle_index = p;
                        const Particle guess = trial[p];
                        trial[p] = algorithmsECCAPFIM::advance_particle(
                            initial[p], electric, magnetic, kinetic_thread_current_[omp_get_thread_num()],
                            sp.is_neutral() ? 0 : sp.charge() / sp.mass(), sp.is_neutral() ? 0 : sp.charge() * sp.mpw(),
                            dt_, domain, particle_tolerance_, particle_iterations_, &bc_handler, policy,
                            &state.events(cell)[p], state.warm ? &guess : nullptr);
                    }
                } catch (...) {
#pragma omp critical(eccapfim_particle_failure)
                    {
                        if (!failure) {
                            failure = std::current_exception();
                            failed_cell = cell;
                            failed_particle = particle_index;
                        }
                    }
                }
            }
        }
        if (failure) {
            // Only the simulation thread writes files, after all particle workers join.
            nlohmann::json report{{"step", current_timestep_},
                                  {"field_iteration", nonlinear_iterations},
                                  {"species", name},
                                  {"cell", failed_cell},
                                  {"particle_index", failed_particle},
                                  {"dt", dt_},
                                  {"system_config", system_config},
                                  {"particles_config", particles_config},
                                  {"field_format", "3 int32 dimensions followed by doubles; Field3d storage order"}};
            const auto& particle = state.initial(failed_cell)[failed_particle];
            report["initial"] = {{"position", {particle.coord.x(), particle.coord.y(), particle.coord.z()}},
                                 {"velocity", {particle.velocity.x(), particle.velocity.y(), particle.velocity.z()}}};
            try {
                std::rethrow_exception(failure);
            } catch (const algorithmsECCAPFIM::ParticleSolveFailure& error) {
                report["particle_solver"] = error.details;
                report["error"] = error.what();
            } catch (const std::exception& error) {
                report["error"] = error.what();
            }
            std::cerr << "ECCAPFIM_FAILURE " << report.dump() << std::endl;
            if (failure_dump_) {
                try {
                    const auto directory =
                        std::filesystem::path("EccapfimFailure") / ("step_" + std::to_string(current_timestep_));
                    std::filesystem::create_directories(directory);
                    write_field_to_file((directory / "electric.backup").string(), electric);
                    write_field_to_file((directory / "magnetic.backup").string(), magnetic);
                    std::ofstream file;
                    file.exceptions(std::ios::failbit | std::ios::badbit);
                    file.open(directory / "failure.json");
                    file << report.dump(2) << '\n';
                    file.close();
                    std::cerr << "ECCAPFIM failure dump: " << directory << std::endl;
                } catch (const std::exception& error) {
                    std::cerr << "ECCAPFIM failure dump write failed: " << error.what() << std::endl;
                }
            }
            std::rethrow_exception(failure);
        }
        state.warm = true;
        {
            RECORD_TIMER_NAMED("eccapfim current fold");
#pragma omp parallel for schedule(static)
            for (int p = 0; p < state.current.capacity(); ++p) {
                double sum = 0.0;
                for (int t = 0; t < threads; ++t) sum += kinetic_thread_current_[t](p);
                state.current(p) = sum;
            }
            algorithmsECCAPFIM::copy_periodic_nodes(state.current, domain, &bc_handler);
            fieldJ += state.current;
        }
    }
}

// `unknown_range_` and the arguments are storage indices (with ghost layers),
// the same numbering as pack().
int SimulationEccapfim::reduced_dof(int i, int j, int k, int component) const {
    const int storage[3] = {i, j, k};
    int index[3];
    for (int d = 0; d < 3; ++d) {
        index[d] = bc_handler.is_periodic(d) ? domain.grid.periodic_node_owner(storage[d], d) : storage[d];
        if (index[d] < unknown_range_.start[d] || index[d] >= unknown_range_.end[d])
            return -1;
    }
    const auto n = unknown_range_.end - unknown_range_.start;
    return 3 * ((index[0] - unknown_range_.start.x()) * n.y() * n.z() + (index[1] - unknown_range_.start.y()) * n.z() +
                (index[2] - unknown_range_.start.z())) +
           component;
}

void SimulationEccapfim::begin_nonlinear_step() {
    kinetic_states_ = save_particles();
}

void SimulationEccapfim::evaluate_current(const Field3d& electric, const Field3d& magnetic) {
    push_particles_and_deposit_current(kinetic_states_, electric, magnetic);
}

void SimulationEccapfim::accept_nonlinear_step() {
    accept_particles(kinetic_states_);
}

void SimulationEccapfim::solve_electric_field() {
    RECORD_TIMER;
    const Field3d magnetic = fieldB + fieldBInit;   // Lorentz force uses B^n + B_external.
    const Vector rhs = pack(fieldE + (0.5 * dt_) * (maxwell.curlB * fieldB)).cwiseProduct(electric_free_);
    auto residual = [&](const Vector& electric) -> Vector {
        if (!electric.allFinite())
            throw std::runtime_error("ECCAPFIM: nonfinite field iterate");
        unpack(electric, fieldEp);
        const Field3d full_electric = fieldEp + fieldE_external;
        ++residual_evaluations_;
        evaluate_current(full_electric, magnetic);
        // Maxwell residual F(Eh) = M Eh - E^n - dt/2 curl B^n + dt/2 J[Eh].
        return maxwell_ * electric - rhs + (0.5 * dt_) * pack(fieldJ).cwiseProduct(electric_free_);
    };

    Vector x = pack(fieldE);
    std::vector<Vector> old_f, old_g;
    std::vector<Vector> old_x, old_r;
    double initial_norm = 0;
    bool converged = false;
    for (int it = 0; it <= max_iterations_; ++it) {
        nonlinear_iterations = it;
        Vector r;
        {
            RECORD_TIMER_NAMED("eccapfim residual");
            r = residual(x);
        }
        nonlinear_residual = r.norm();
        if (iteration_log_)
            std::cout << "ECCAPFIM_RESIDUAL step=" << current_timestep_ << " iteration=" << it
                      << " norm=" << std::setprecision(17) << nonlinear_residual << " target="
                      << absolute_tolerance_ + relative_tolerance_ * (it == 0 ? nonlinear_residual : initial_norm)
                      << std::endl;
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
        Vector f;
        bool maxwell_ok = false;
        size_t maxwell_iterations = 0;
        {
            RECORD_TIMER_NAMED("eccapfim Maxwell solve");
            if (use_amgcl_) {
                Vector solution = Vector::Zero(r.size());
                std::tie(amgcl_iterations_, amgcl_error_) = (*amgcl_solver_)(r, solution);
                f = mixing_ * solution;
                maxwell_iterations = amgcl_iterations_;
                maxwell_ok = std::isfinite(amgcl_error_) && amgcl_iterations_ < 2000;
            } else if (use_ic_) {
                f = mixing_ * maxwell_ic_solver_.solve(r);
                maxwell_iterations = maxwell_ic_solver_.iterations();
                maxwell_ok = maxwell_ic_solver_.info() == Eigen::Success;
            } else if (use_mass_matrix_) {
                if (mass_matrix_every_iteration_)
                    mass_matrix_ready_ = false;
                ensure_mass_matrix(magnetic);
                f = mixing_ * solve_mass_system(r);
                maxwell_iterations = 0;
                maxwell_ok = true;
            } else if (use_block_jacobi_) {
                f = mixing_ * maxwell_block_solver_.solve(r);
                maxwell_iterations = maxwell_block_solver_.iterations();
                maxwell_ok = maxwell_block_solver_.info() == Eigen::Success;
            } else {
                f = mixing_ * maxwell_solver_.solve(r);
                maxwell_iterations = maxwell_solver_.iterations();
                maxwell_ok = maxwell_solver_.info() == Eigen::Success;
            }
        }
        if (!maxwell_ok || !f.allFinite())
            throw std::runtime_error("ECCAPFIM: Maxwell preconditioner did not converge");
        if (iteration_log_)
            std::cout << "ECCAPFIM_MAXWELL step=" << current_timestep_ << " iteration=" << it
                      << " iterations=" << maxwell_iterations << std::endl;
        const Vector g = x - f;
        Vector next = g;
        if (use_ngmres_) {
            // Nonlinear GMRES on the right-preconditioned residual history:
            // x <- x + dX c, where dX and dR are differences of previous
            // accepted iterates and residuals.
            const int m = static_cast<int>(old_x.size());
            if (m > 0) {
                Eigen::MatrixXd dX(x.size(), m), dR(x.size(), m);
                for (int j = 0; j < m; ++j) {
                    dX.col(j) = old_x[j] - x;
                    dR.col(j) = old_r[j] - r;
                }
                const Vector coefficients = dR.colPivHouseholderQr().solve(-r);
                const Vector accelerated = g + dX * coefficients;
                if (accelerated.allFinite() && (accelerated - g).norm() <= 10 * f.norm())
                    next = accelerated;
                else {
                    old_x.clear();
                    old_r.clear();
                }
            }
            if (history_ > 0) {
                old_x.push_back(x);
                old_r.push_back(r);
                if (static_cast<int>(old_x.size()) > history_) {
                    old_x.erase(old_x.begin());
                    old_r.erase(old_r.begin());
                }
            }
        } else {
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
        }
        x = next;
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: nonlinear field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
}


void SimulationEccapfim::accept_particles(ParticleIterations& states) {
    RECORD_TIMER;
    for (auto& [name, state] : states) {
        auto& sp = *state.species;
        // Full rebin: implicit orbits may cross several cells or the periodic seam.
        for (int cell = 0; cell < sp.size(); ++cell) sp(cell).clear();
        for (int cell = 0; cell < state.trial.capacity(); ++cell)
            for (size_t p = 0; p < state.trial(cell).size(); ++p) {
                const auto& particle = state.trial(cell)[p];
                const auto& event = state.events(cell)[p];
                if (event.removed != Face::UNAVAILABLE)
                    sp.step_state.add_loss(event.removed, get_energy_particle(particle.velocity, sp.mass(), sp.mpw()));
                else
                    sp.add_particle(particle);
                for (size_t f = 0; f < event.reflected.size(); ++f)
                    if (event.reflected[f] != 0)
                        sp.step_state.boundary[static_cast<Face>(f)].reflected_count += event.reflected[f];
            }
        sp.currentOnGrid = state.current;
    }
}

void SimulationEccapfim::make_step(int timestep) {
    current_timestep_ = timestep;
    RECORD_TIMER;
    timer::timer timerTotal("Total");
    // x^n and v^n remain unchanged while Eh, particle orbits and J are iterated.
    if (iteration_log_)
        std::cout << "ECCAPFIM_BEGIN step=" << timestep << std::endl;
    // The mass matrix belongs to this step: particles and B have moved.
    mass_matrix_ready_ = false;
    residual_evaluations_ = 0;
    mass_solve_iterations_ = 0;
    begin_nonlinear_step();
    if (coupled_fields_) {
        solve_fields_coupled();
        accept_nonlinear_step();
        timerTotal.finish();
        std::cout << "ECCAPFIM coupled iterations=" << nonlinear_iterations
                  << " residual=" << nonlinear_residual << '\n';
        return;
    }
    solve_electric_field();

    // Eh = (E^{n+1} + E^n)/2; Faraday advances the self-consistent B field.
    fieldEn = 2.0 * fieldEp - fieldE;
    maxwell.advance_magnetic_field(fieldBn, fieldB, fieldE, fieldEn, dt_);
    algorithmsECCAPFIM::copy_periodic_nodes(fieldEn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    algorithmsECCAPFIM::copy_periodic_nodes(fieldBn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);

    // The last residual used the accepted Eh: commit those orbits and their J once.
    accept_nonlinear_step();
    timerTotal.finish();
    std::cout << "ECCAPFIM iterations=" << nonlinear_iterations << " residual=" << nonlinear_residual << '\n';
}

Vector3R SimulationEccapfim::field_momentum(const Field3d& electric, const Field3d& magnetic) const {
    // Cell (i, j, k) has its centre at (i + 1/2, j + 1/2, k + 1/2): E_x at
    // (i + 1/2, j, k) is averaged over the four edges, B_x at (i, j + 1/2, k + 1/2)
    // over the two faces, and cyclically for y and z.
    const Vector3I n = domain.num_cells();
    double px = 0, py = 0, pz = 0;
#pragma omp parallel for collapse(2) reduction(+ : px, py, pz)
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k) {
                int index[3][2];
                const int cell[3] = {i, j, k};
                for (int a = 0; a < 3; ++a)
                    for (int s = 0; s < 2; ++s)
                        index[a][s] = algorithmsECCAPFIM::node_owner(cell[a] + s, a, domain, &bc_handler);
                Vector3R E(0, 0, 0), B(0, 0, 0);
                for (int c = 0; c < 3; ++c) {
                    const int b = (c + 1) % 3, d = (c + 2) % 3;
                    for (int s = 0; s < 2; ++s)
                        for (int t = 0; t < 2; ++t) {
                            int e[3];
                            e[c] = index[c][0];
                            e[b] = index[b][s];
                            e[d] = index[d][t];
                            E[c] += 0.25 * electric(e[0], e[1], e[2], c);
                        }
                    for (int s = 0; s < 2; ++s) {
                        int m[3];
                        m[c] = index[c][s];
                        m[b] = index[b][0];
                        m[d] = index[d][0];
                        B[c] += 0.5 * magnetic(m[0], m[1], m[2], c);
                    }
                }
                const Vector3R density = E.cross(B);
                px += density.x();
                py += density.y();
                pz += density.z();
            }
    return Vector3R(px, py, pz);
}

void SimulationEccapfim::make_diagnostic(int timestep) {
    RECORD_TIMER;
    if (!diagnostic_) {
        diagnostic_ = std::make_unique<Diagnostics>(species);
        std::vector<std::string> guiding_centers;
        if (const MomentSources* extra = extra_moment_sources())
            for (const auto& entry : *extra) guiding_centers.push_back(entry.first);
        outputs_ = OutputFactory::create(system_config.value("diagnostics", nlohmann::json::object()), domain, species,
                                         system_config, fieldEn, fieldBn, fieldE_external, fieldBInit, guiding_centers);
    }
    auto& diag = *diagnostic_;
    const auto range = bc_handler.active_range(domain.grid);
    double kinetic = 0, previous = 0, injected = 0, external_work = 0, lost = 0;
    for (const auto& [name, ptr] : species) {
        const auto moments = ptr->diagnostic_moments();
        diag.addEnergy(name, moments.energy());
        diag.addEnergy(name + "Particles", static_cast<double>(moments.count));
        diag.addEnergy(name + "Inject", ptr->step_state.injection_energy);
        diag.addEnergy(name + "XY", moments.kinetic.x() + moments.kinetic.y());
        diag.addEnergy(name + "Z", moments.kinetic.z());
        diag.addBoundary(name + "_InjectE", ptr->step_state.injection_energy);
        diag.addBoundary(name + "_InjectN", ptr->step_state.injection_count);
        const char* face_names[] = {"XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"};
        int fi = 0;
        for (const auto face : {Face::XMIN, Face::XMAX, Face::YMIN, Face::YMAX, Face::ZMIN, Face::ZMAX}) {
            const auto it = ptr->step_state.boundary.find(face);
            const BoundaryStepEvents stats = it != ptr->step_state.boundary.end() ? it->second : BoundaryStepEvents{};
            const std::string prefix = name + "_" + face_names[fi++];
            diag.addBoundary(prefix + "LostE", stats.lost_energy);
            diag.addBoundary(prefix + "LostN", stats.lost_count);
            diag.addBoundary(prefix + "ReflN", stats.reflected_count);
            diag.addBoundary(prefix + "EmitE", 0);
            diag.addBoundary(prefix + "EmitN", 0);
        }
        lost += ptr->step_state.total_lost_energy();
        kinetic += moments.energy();
        previous += ptr->step_state.kinetic_energy_before_step;
        injected += ptr->step_state.injection_energy;
        // currentOnGrid is produced by the first accepted step; at t=0 it is
        // not yet a physical current and must not contribute to diagnostics.
        if (timestep > 0)
            external_work += dt_ * dot_product_sum(fieldE_external, ptr->currentOnGrid, range);
    }
    kinetic += extra_kinetic_energy();
    previous += extra_previous_kinetic_energy();
    if (timestep > 0)
        external_work += extra_external_work();
    add_extra_diagnostics(diag);
    const double e = (0.5 * dot_product_sum(fieldEn, fieldEn, range)), b = (0.5 * dot_product_sum(fieldBn, fieldBn, range));
    const double delta = e + b - (0.5 * dot_product_sum(fieldE, fieldE, range)) - (0.5 * dot_product_sum(fieldB, fieldB, range));
    // Magnetization enters through J = J_free + curl_B M, so its work M . dB is
    // already part of E . J and of the field energy change.
    const double balance = kinetic - previous - injected + delta - external_work + lost;
    diag.addEnergy("energyFieldE", e);
    diag.addEnergy("energyFieldB", b);
    diag.addEnergy("energyFieldBFull", (0.5 * dot_product_sum(fieldBn + fieldBInit, fieldBn + fieldBInit, range)));
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
    diag.addEnergy("nonlinearEvaluations", residual_evaluations_);
    diag.addEnergy("massSolveIterations", static_cast<double>(mass_solve_iterations_));
    // Momentum of the accepted step: sum m w v of the markers (guiding centres:
    // m w v_par b) and sum E x B of the self fields. The scheme conserves
    // energy, not momentum: the columns show the momentum error.
    Vector3R particle_momentum = extra_momentum();
    for (const auto& [name, ptr] : species)
        for (int cell = 0; cell < ptr->size(); ++cell)
            for (const Particle& p : (*ptr)(cell)) particle_momentum += (ptr->mass() * ptr->mpw()) * p.velocity;
    const Vector3R momentum_field = field_momentum(fieldEn, fieldBn);
    const char* axis_name[3] = {"X", "Y", "Z"};
    for (int d = 0; d < 3; ++d) {
        diag.addEnergy(std::string("momentumParticles") + axis_name[d], particle_momentum[d]);
        diag.addEnergy(std::string("momentumField") + axis_name[d], momentum_field[d]);
        diag.addEnergy(std::string("momentumTotal") + axis_name[d], particle_momentum[d] + momentum_field[d]);
    }
    write_diagnostics(timestep, diag,
                      {domain, species, fieldEn, fieldBn, fieldE_external, fieldBInit, bc_handler,
                       diagnostic_density_shape(), extra_moment_sources()},
                      outputs_);
    if (failure_dump_) {
        for (const auto& output :
             system_config.value("diagnostics", nlohmann::json::object()).value("outputs", nlohmann::json::array())) {
            const int interval = output.value("interval", 0);
            if (output.value("type", std::string()) == "recovery" && interval > 0 && timestep % interval == 0) {
                nlohmann::json checkpoint{{"step", timestep},
                                          {"time", timestep * dt_},
                                          {"dt", dt_},
                                          {"phase", "accepted state after collisions"},
                                          {"rng_state_saved", false}};
                for (const auto& [name, sp] : species)
                    checkpoint["particle_counts"][name] = sp->get_total_num_of_particles();
                std::ofstream file;
                file.exceptions(std::ios::failbit | std::ios::badbit);
                file.open("Recovery/checkpoint.json");
                file << checkpoint.dump(2) << '\n';
            }
        }
    }
}

void SimulationEccapfim::finalize() {
    RECORD_TIMER;
    for (auto& output : outputs_) output->finalize();
    Simulation::finalize();
}
