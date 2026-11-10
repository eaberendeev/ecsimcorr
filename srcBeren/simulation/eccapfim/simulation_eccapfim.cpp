#include "simulation_eccapfim.h"

#include <Eigen/Core>
#include <algorithm>
#include <Eigen/QR>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
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
    for (auto& [name, state] : states) {
        auto& sp = *state.species;
        state.current.setZero();
        const auto policy = bc_handler.path_policy(sp, domain);
        std::exception_ptr failure;
        int failed_cell = -1;
        size_t failed_particle = 0;
        {
            RECORD_TIMER_NAMED("eccapfim particle push");
#pragma omp parallel for schedule(dynamic, 64)
            for (int cell = 0; cell < state.initial.capacity(); ++cell) {
                size_t particle_index = 0;
                try {
                    const auto& initial = state.initial(cell);
                    auto& trial = state.trial(cell);
                    for (size_t p = 0; p < initial.size(); ++p) {
                        particle_index = p;
                        const Particle guess = trial[p];
                        trial[p] = algorithmsECCAPFIM::advance_particle(
                            initial[p], electric, magnetic, state.current,
                            sp.is_neutral() ? 0 : sp.charge() / sp.mass(),
                            sp.is_neutral() ? 0 : sp.charge() * sp.mpw(), dt_, domain, particle_tolerance_,
                            particle_iterations_, &bc_handler, policy, &state.events(cell)[p],
                            state.warm ? &guess : nullptr);
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

void SimulationEccapfim::add_particle_mass(const Vector3R& position, const Eigen::Matrix3d& response) {
    // E_c at (i + 1/2, j, k) for c = x, and cyclic: linear weights along c on
    // half-integer nodes and quadratic weights across c on integer nodes.
    const Vector3R r = domain.to_cell_coordinates(position);
    const auto spline1 = [](double x) { return std::max(0.0, 1.0 - std::abs(x)); };
    const auto spline2 = [](double x) {
        x = std::abs(x);
        return x <= 0.5 ? 0.75 - x * x : (x < 1.5 ? 0.5 * (1.5 - x) * (1.5 - x) : 0.0);
    };
    struct Weight {
        int node[3];
        int row;
        double value;
    };
    Weight weights[3][18];
    int count[3] = {0, 0, 0};
    for (int c = 0; c < 3; ++c) {
        int first[3], width[3];
        double w[3][3];
        for (int a = 0; a < 3; ++a) {
            if (a == c) {
                first[a] = static_cast<int>(std::floor(r[a] - 0.5));
                width[a] = 2;
                for (int n = 0; n < 2; ++n) w[a][n] = spline1(r[a] - (first[a] + n + 0.5));
            } else {
                first[a] = static_cast<int>(std::floor(r[a] + 0.5)) - 1;
                width[a] = 3;
                for (int n = 0; n < 3; ++n) w[a][n] = spline2(r[a] - (first[a] + n));
            }
        }
        for (int i = 0; i < width[0]; ++i)
            for (int j = 0; j < width[1]; ++j)
                for (int k = 0; k < width[2]; ++k) {
                    const double value = w[0][i] * w[1][j] * w[2][k];
                    if (value == 0.0)
                        continue;
                    const int node[3] = {first[0] + i, first[1] + j, first[2] + k};
                    const int row = reduced_dof(node[0] + GHOST_CELLS, node[1] + GHOST_CELLS, node[2] + GHOST_CELLS, c);
                    weights[c][count[c]++] = {{node[0], node[1], node[2]}, row, value};
                }
    }
    for (int c1 = 0; c1 < 3; ++c1)
        for (int a = 0; a < count[c1]; ++a) {
            const Weight& wa = weights[c1][a];
            // Constrained rows keep the identity of the Maxwell operator.
            if (wa.row < 0 || electric_free_[wa.row] == 0)
                continue;
            for (int c2 = 0; c2 < 3; ++c2) {
                const double t = response(c1, c2) * wa.value;
                if (t == 0.0)
                    continue;
                for (int b = 0; b < count[c2]; ++b) {
                    const Weight& wb = weights[c2][b];
                    stencil_operator_.add(wa.row, c1, c2, wb.node[0] - wa.node[0], wb.node[1] - wa.node[1],
                                          wb.node[2] - wa.node[2], t * wb.value);
                }
            }
        }
}

void SimulationEccapfim::assemble_mass_matrix(const Field3d& magnetic) {
    // Crank-Nicolson response: v_h = R(beta) (v^n + alpha E), alpha = dt q / 2m,
    // beta = alpha B, so (dt/2) dJ/dE = (dt/2)^2 q^2 w / m R(beta).
    for (const auto& [name, state] : kinetic_states_) {
        const auto& sp = *state.species;
        if (sp.is_neutral())
            continue;
        const double q = sp.charge(), m = sp.mass(), mpw = sp.mpw();
#pragma omp parallel for schedule(dynamic, 16)
        for (int cell = 0; cell < state.trial.capacity(); ++cell)
            for (size_t p = 0; p < state.trial(cell).size(); ++p) {
                // The orbit kernels average E along x^n -> x^{n+1}; its midpoint is the
                // closest point kernel. Orbits split by wall events are rare and use it too.
                const Vector3R midpoint = mass_matrix_point(state.initial(cell)[p].coord, state.trial(cell)[p].coord);
                const Vector3R B = algorithmsECCAPFIM::gather_magnetic_point(magnetic, midpoint, domain, &bc_handler);
                const Vector3R beta = B * (0.5 * dt_ * q / m);
                const double denom = 1.0 + beta.squaredNorm();
                const double a = 0.25 * dt_ * dt_ * q * q * mpw / m / denom;
                Eigen::Matrix3d T;
                T << a * (1 + beta.x() * beta.x()), a * (beta.z() + beta.x() * beta.y()),
                    a * (-beta.y() + beta.x() * beta.z()), a * (-beta.z() + beta.y() * beta.x()),
                    a * (1 + beta.y() * beta.y()), a * (beta.x() + beta.y() * beta.z()),
                    a * (beta.y() + beta.z() * beta.x()), a * (-beta.x() + beta.z() * beta.y()),
                    a * (1 + beta.z() * beta.z());
                add_particle_mass(algorithmsECCAPFIM::wrap_position(midpoint, domain, &bc_handler), T);
            }
    }
}

Vector3R SimulationEccapfim::mass_matrix_point(const Vector3R& start, const Vector3R& end) const {
    if (mass_matrix_at_end_)
        return end;
    // The trial end point may be wrapped through a periodic face: undo that shift.
    Vector3R displacement = end - start;
    for (int d = 0; d < 3; ++d)
        if (bc_handler.is_periodic(d)) {
            const double length = domain.num_cells(d) * domain.cell_size(d);
            displacement[d] -= length * std::round(displacement[d] / length);
        }
    return start + 0.5 * displacement;
}

// A = maxwell_ + (dt/2) M_p in the fixed particle pattern; GMRES with node
// block-Jacobi by default (xpic: PETSc GMRES, rtol 1e-6), or sparse LU.
void SimulationEccapfim::rebuild_mass_matrix(const Field3d& magnetic) {
    RECORD_TIMER;
    {
        RECORD_TIMER_NAMED("mass matrix assembly");
        stencil_operator_.reset();
        assemble_mass_matrix(magnetic);
    }
    timer::timer factorization_timer("mass matrix factorization");
    if (mass_matrix_lu_) {
        mass_preconditioner_ = stencil_operator_.to_sparse();
        mass_solver_.analyzePattern(mass_preconditioner_);
        mass_solver_.factorize(mass_preconditioner_);
        if (mass_solver_.info() != Eigen::Success)
            throw std::runtime_error("ECCAPFIM: mass-matrix factorization failed: " + mass_solver_.lastErrorMessage());
    } else {
        stencil_operator_.factor(mass_preconditioner_kind_);
    }
    mass_matrix_ready_ = true;
}

// Built once per step at the first use; the solvers clear the flag at every
// outer iteration when `mass_matrix_every_iteration` is set.
void SimulationEccapfim::ensure_mass_matrix(const Field3d& magnetic) {
    if (!mass_matrix_ready_)
        rebuild_mass_matrix(magnetic);
}

SimulationEccapfim::Vector SimulationEccapfim::solve_mass_system(const Vector& rhs) {
    RECORD_TIMER;
    Vector solution;
    if (mass_matrix_lu_) {
        solution = mass_solver_.solve(rhs);
        if (mass_solver_.info() != Eigen::Success)
            throw std::runtime_error("ECCAPFIM: mass-matrix solve failed");
    } else {
        int iterations = 0;
        if (!stencil_operator_.solve_gmres(rhs, solution, mass_matrix_tolerance_, 50, 2000, iterations))
            throw std::runtime_error("ECCAPFIM: mass-matrix GMRES did not converge");
        mass_solve_iterations_ += iterations;
    }
    if (!solution.allFinite())
        throw std::runtime_error("ECCAPFIM: nonfinite mass-matrix solution");
    return solution;
}

Eigen::SparseMatrix<double> SimulationEccapfim::mass_only_matrix() const {
    return stencil_operator_.to_sparse() - maxwell_;
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


SimulationEccapfim::Vector SimulationEccapfim::pack6(const Field3d& electric, const Field3d& magnetic) const {
    const auto n = unknown_range_.end - unknown_range_.start;
    const auto first = unknown_range_.start;
    Vector result(6 * n.elements_product());
    int p = 0;
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k) {
                for (int c = 0; c < 3; ++c) result[p++] = electric(i + first.x(), j + first.y(), k + first.z(), c);
                for (int c = 0; c < 3; ++c) result[p++] = magnetic(i + first.x(), j + first.y(), k + first.z(), c);
            }
    return result;
}

void SimulationEccapfim::unpack6(const Vector& vector, Field3d& electric, Field3d& magnetic) const {
    const auto n = unknown_range_.end - unknown_range_.start;
    const auto first = unknown_range_.start;
    int p = 0;
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k) {
                for (int c = 0; c < 3; ++c) electric(i + first.x(), j + first.y(), k + first.z(), c) = vector[p++];
                for (int c = 0; c < 3; ++c) magnetic(i + first.x(), j + first.y(), k + first.z(), c) = vector[p++];
            }
    algorithmsECCAPFIM::copy_periodic_nodes(electric, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(magnetic, domain, &bc_handler);
    bc_handler.apply_to_fields(electric, FieldType::ELECTRIC, domain);
    bc_handler.apply_to_fields(magnetic, FieldType::MAGNETIC, domain);
}

SimulationEccapfim::Vector SimulationEccapfim::coupled_residual(const Vector& state) {
    if (!state.allFinite())
        throw std::runtime_error("ECCAPFIM: nonfinite coupled field iterate");
    unpack6(state, fieldEp, fieldBHalf);
    const Field3d full_electric = fieldEp + fieldE_external;
    const Field3d full_magnetic = fieldBHalf + fieldBInit;
    ++residual_evaluations_;
    evaluate_current(full_electric, full_magnetic);

    // F_E = E_h - E^n + dt/2 (J - curl_B B_h),
    // F_B = B_h - B^n + dt/2 curl_E E_h.
    // `fieldJ` already contains the kinetic current and curl(M) for the
    // drift-kinetic backend.
    Field3d residual_e = fieldEp - fieldE + (0.5 * dt_) * (fieldJ - (maxwell.curlB * fieldBHalf));
    Field3d residual_b = fieldBHalf - fieldB + (0.5 * dt_) * (maxwell.curlE * fieldEp);
    algorithmsECCAPFIM::copy_periodic_nodes(residual_e, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(residual_b, domain, &bc_handler);
    bc_handler.apply_to_fields(residual_e, FieldType::ELECTRIC, domain);
    bc_handler.apply_to_fields(residual_b, FieldType::MAGNETIC, domain);
    return pack6(residual_e, residual_b);
}

SimulationEccapfim::Vector SimulationEccapfim::jacobian_vector(const Vector& state, const Vector& direction,
                                                              const Vector& base_residual) {
    const double direction_norm = direction.norm();
    if (direction_norm == 0.0)
        return Vector::Zero(base_residual.size());
    const double epsilon = 1e-7 * (1.0 + state.norm()) / direction_norm;
    const Vector perturbed = state + epsilon * direction;
    return (coupled_residual(perturbed) - base_residual) / epsilon;
}

SimulationEccapfim::Vector SimulationEccapfim::apply_right_preconditioner(const Vector& vector) {
    if (!use_mass_matrix_)
        return vector;
    // fieldBHalf holds the last evaluated iterate: B^n at the first call.
    ensure_mass_matrix(fieldBHalf + fieldBInit);
    Field3d field_e(domain.size(), 3), field_b(domain.size(), 3);
    unpack6(vector, field_e, field_b);
    const Field3d rhs_e = field_e + (0.5 * dt_) * (maxwell.curlB * field_b);
    const Vector rhs_vector = pack(rhs_e);
    const Vector correction_e = solve_mass_system(rhs_vector);
    Field3d correction_e_field(domain.size(), 3);
    unpack(correction_e, correction_e_field);
    Field3d correction_b = field_b - (0.5 * dt_) * (maxwell.curlE * correction_e_field);
    algorithmsECCAPFIM::copy_periodic_nodes(correction_e_field, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(correction_b, domain, &bc_handler);
    bc_handler.apply_to_fields(correction_e_field, FieldType::ELECTRIC, domain);
    bc_handler.apply_to_fields(correction_b, FieldType::MAGNETIC, domain);
    return pack6(correction_e_field, correction_b);
}

SimulationEccapfim::Vector SimulationEccapfim::gmres(const Vector& state, const Vector& rhs, int& iterations,
                                                     bool& ok) {
    const int n = rhs.size();
    iterations = 0;
    ok = true;
    if (n == 0)
        return Vector::Zero(0);
    const double beta = rhs.norm();
    if (beta == 0.0)
        return Vector::Zero(n);

    const int max_iterations = std::min(coupled_gmres_max_, n);
    std::vector<Vector> basis(max_iterations + 1);
    Eigen::MatrixXd hessenberg = Eigen::MatrixXd::Zero(max_iterations + 1, max_iterations);
    basis[0] = rhs / beta;
    const Vector base_residual = -rhs;   // rhs = -F(state).
    Vector solution = Vector::Zero(n);
    double linear_residual = beta;

    for (int j = 0; j < max_iterations; ++j) {
        const Vector preconditioned = apply_right_preconditioner(basis[j]);
        Vector w = jacobian_vector(state, preconditioned, base_residual);
        for (int i = 0; i <= j; ++i) {
            hessenberg(i, j) = w.dot(basis[i]);
            w -= hessenberg(i, j) * basis[i];
        }
        // One reorthogonalization pass keeps the small Krylov basis stable.
        for (int i = 0; i <= j; ++i) {
            const double correction = w.dot(basis[i]);
            hessenberg(i, j) += correction;
            w -= correction * basis[i];
        }
        hessenberg(j + 1, j) = w.norm();
        if (hessenberg(j + 1, j) > 1e-15)
            basis[j + 1] = w / hessenberg(j + 1, j);

        Eigen::VectorXd rhs_small = Eigen::VectorXd::Zero(j + 2);
        rhs_small[0] = beta;
        const Eigen::VectorXd coefficients =
            hessenberg.topLeftCorner(j + 2, j + 1).colPivHouseholderQr().solve(rhs_small);
        solution.setZero();
        for (int i = 0; i <= j; ++i) solution += coefficients[i] * basis[i];
        linear_residual = (rhs_small - hessenberg.topLeftCorner(j + 2, j + 1) * coefficients).norm();
        iterations = j + 1;
        if (!solution.allFinite()) {
            ok = false;
            break;
        }
        if (linear_residual <= 1e-10 * beta)
            break;
        if (hessenberg(j + 1, j) <= 1e-15)
            break;
    }
    ok = ok && solution.allFinite();
    return apply_right_preconditioner(solution);
}

// Newton-Krylov: GMRES on finite-difference Jacobian products, each of which
// pushes all particles, with the mass matrix as right linear preconditioner.
void SimulationEccapfim::solve_fields_newton_krylov(Vector& state, Vector& residual) {
    const double initial_norm = residual.norm();
    bool converged = false;
    for (int outer = 0; outer <= coupled_max_outer_; ++outer) {
        nonlinear_iterations = outer;
        nonlinear_residual = residual.norm();
        if (iteration_log_)
            std::cout << "ECCAPFIM_COUPLED step=" << current_timestep_ << " newton=" << outer
                      << " residual=" << std::setprecision(17) << nonlinear_residual << std::endl;
        if (!residual.allFinite() || !std::isfinite(nonlinear_residual))
            throw std::runtime_error("ECCAPFIM: nonfinite coupled residual");
        if (nonlinear_residual <= coupled_atol_ + coupled_rtol_ * initial_norm) {
            converged = true;
            break;
        }
        if (outer == coupled_max_outer_)
            break;

        if (mass_matrix_every_iteration_)
            mass_matrix_ready_ = false;
        int linear_iterations = 0;
        bool linear_ok = false;
        const Vector direction = gmres(state, -residual, linear_iterations, linear_ok);
        if (!linear_ok || !direction.allFinite())
            throw std::runtime_error("ECCAPFIM: coupled GMRES failed");
        if (iteration_log_)
            std::cout << "ECCAPFIM_COUPLED_GMRES step=" << current_timestep_ << " newton=" << outer
                      << " iterations=" << linear_iterations << std::endl;

        double alpha = 1.0;
        Vector trial_state, trial_residual;
        bool accepted = false;
        for (int line_search = 0; line_search < 10; ++line_search) {
            trial_state = state + alpha * direction;
            trial_residual = coupled_residual(trial_state);
            if (trial_residual.norm() < (1.0 - 1e-4 * alpha) * residual.norm()) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted)
            throw std::runtime_error("ECCAPFIM: coupled Newton line search failed");
        state = trial_state;
        residual = trial_residual;
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: coupled field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
}

// Nonlinear GMRES with right nonlinear preconditioning (Oosterlee & Washio
// 2000; PETSc SNESNGMRES as used by xpic). One preconditioned step
// x^M = x - P^{-1} F(x) uses the frozen-response Jacobian
// P = [[I + dt/2 M, -dt/2 curl_B], [dt/2 curl_E, I]]; then x^A minimizes the
// linearized residual over the history of iterates. Every candidate is judged
// by the full nonlinear residual, i.e. by a complete particle push.
void SimulationEccapfim::solve_fields_ngmres(Vector& state, Vector& residual) {
    constexpr double gamma_a = 2.0, epsilon_b = 0.1, delta_b = 0.9, gamma_c = 2.0;
    constexpr int restart_after = 2;
    const double initial_norm = residual.norm();
    Vector last_evaluated = state;
    std::vector<Vector> history_x, history_f;
    std::vector<double> history_norm;
    int restart_count = 0;
    bool converged = false;
    for (int outer = 0; outer <= coupled_max_outer_; ++outer) {
        nonlinear_iterations = outer;
        nonlinear_residual = residual.norm();
        if (iteration_log_)
            std::cout << "ECCAPFIM_NGMRES step=" << current_timestep_ << " iteration=" << outer
                      << " residual=" << std::setprecision(17) << nonlinear_residual << std::endl;
        if (!residual.allFinite() || !std::isfinite(nonlinear_residual))
            throw std::runtime_error("ECCAPFIM: nonfinite coupled residual");
        if (nonlinear_residual <= coupled_atol_ + coupled_rtol_ * initial_norm) {
            converged = true;
            break;
        }
        if (outer == coupled_max_outer_)
            break;

        history_x.push_back(state);
        history_f.push_back(residual);
        history_norm.push_back(nonlinear_residual);
        if (static_cast<int>(history_x.size()) > coupled_history_) {
            history_x.erase(history_x.begin());
            history_f.erase(history_f.begin());
            history_norm.erase(history_norm.begin());
        }

        // Preconditioned step. fieldBHalf holds B of `state` (last evaluation).
        if (mass_matrix_every_iteration_)
            mass_matrix_ready_ = false;
        if (mass_matrix_check_ && outer == 0)
            check_mass_matrix(state, residual);
        const Vector state_m = state - apply_right_preconditioner(residual);
        const Vector residual_m = coupled_residual(state_m);
        last_evaluated = state_m;
        const double norm_m = residual_m.norm();

        // Least-squares combination over the history, linearized around x^M.
        const int m = static_cast<int>(history_x.size());
        Eigen::MatrixXd dF(residual.size(), m);
        for (int j = 0; j < m; ++j) dF.col(j) = history_f[j] - residual_m;
        const Vector beta = dF.colPivHouseholderQr().solve(-residual_m);
        Vector state_a = state_m;
        for (int j = 0; j < m; ++j) state_a += beta[j] * (history_x[j] - state_m);
        bool select_a = false;
        Vector residual_a;
        double norm_a = std::numeric_limits<double>::infinity();
        const double min_norm = *std::min_element(history_norm.begin(), history_norm.end());
        if (state_a.allFinite() && (state_a - state_m).norm() > 0.0) {
            residual_a = coupled_residual(state_a);
            last_evaluated = state_a;
            norm_a = residual_a.norm();
            if (norm_a * norm_a <= gamma_a * min_norm * min_norm) {
                select_a = true;
            } else {
                double min_distance = std::numeric_limits<double>::infinity();
                for (const auto& x : history_x) min_distance = std::min(min_distance, (state_a - x).norm());
                select_a = epsilon_b * (state_a - state_m).norm() < min_distance ||
                           std::sqrt(norm_a) < delta_b * std::sqrt(min_norm);
            }
            select_a = select_a && std::isfinite(norm_a) && norm_a <= norm_m * 10;
        }
        // Restart the history after repeated growth of the accelerated residual.
        restart_count = (std::isfinite(norm_a) && norm_a > gamma_c * min_norm) ? restart_count + 1 : 0;
        if (restart_count >= restart_after) {
            history_x.clear();
            history_f.clear();
            history_norm.clear();
            restart_count = 0;
        }
        if (select_a) {
            state = state_a;
            residual = residual_a;
        } else {
            state = state_m;
            residual = residual_m;
        }
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: coupled NGMRES field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
    // Particles and sources must belong to the accepted fields.
    if (last_evaluated != state)
        residual = coupled_residual(state);
}

// Diagnostic of the frozen-response approximation (xpic
// -npc_mass_matrix_check): compares P v with the finite-difference Jacobian
// product (F(x + eps v) - F(x)) / eps for random E and B directions.
void SimulationEccapfim::check_mass_matrix(const Vector& state, const Vector& residual) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> u(-1, 1);
    for (int part = 0; part < 2; ++part) {
        Vector v = Vector::Zero(state.size());
        for (int i = 0; i < v.size(); ++i)
            if ((i % 6 < 3) == (part == 0))
                v[i] = u(rng);
        // P v with P = [[maxwell + mass - dt^2/4 curl_B curl_E, -dt/2 curl_B], [dt/2 curl_E, I]].
        ensure_mass_matrix(fieldBHalf + fieldBInit);
        Field3d ve(domain.size(), 3), vb(domain.size(), 3);
        ve.setZero();
        vb.setZero();
        unpack6(v, ve, vb);
        const Eigen::SparseMatrix<double> mass_matrix = mass_only_matrix();
        const Vector pe = pack(ve) + mass_matrix * pack(ve);
        Field3d pe_field(domain.size(), 3);
        pe_field.setZero();
        unpack(pe, pe_field);
        const Field3d p_e = pe_field - (0.5 * dt_) * (maxwell.curlB * vb);
        const Field3d p_b = vb + (0.5 * dt_) * (maxwell.curlE * ve);
        const Vector pv = pack6(p_e, p_b);
        const double eps = 1e-7 * (1.0 + state.norm()) / v.norm();
        const Vector jv = (coupled_residual(state + eps * v) - residual) / eps;
        // The particle part alone: (dt/2) dJ v against the mass matrix.
        Field3d mv_field(domain.size(), 3);
        mv_field.setZero();
        unpack(mass_matrix * pack(ve), mv_field);
        Field3d zero(domain.size(), 3);
        zero.setZero();
        const Vector mv = pack6(mv_field, zero);
        const Vector response = jv - (pv - mv);
        std::cout << "MASS_MATRIX_CHECK step=" << current_timestep_ << " direction=" << (part == 0 ? "E" : "B")
                  << " |Jv|=" << jv.norm() << " |Pv-Jv|/|Jv|=" << (pv - jv).norm() / jv.norm()
                  << " |dt/2 dJ|=" << response.norm() << " |Mv - dt/2 dJ|/|dt/2 dJ|="
                  << (mv - response).norm() / std::max(response.norm(), 1e-300) << std::endl;
    }
    coupled_residual(state);   // restore particles and sources of `state`
}

// Anderson acceleration of the preconditioned fixed point
// G(x) = x - P^{-1} F(x) (Walker & Ni 2011): one particle push per
// iteration, against two for NGMRES. With history m the update is
// x+ = x + g - (dX + dG) gamma, gamma = argmin |g - dG gamma|, g = G(x) - x.
// Convergence and acceptance use the full nonlinear residual F, and the
// accepted state is always the last evaluated one.
void SimulationEccapfim::solve_fields_anderson(Vector& state, Vector& residual) {
    const double initial_norm = residual.norm();
    std::vector<Vector> history_x, history_g;
    bool converged = false;
    for (int outer = 0; outer <= coupled_max_outer_; ++outer) {
        nonlinear_iterations = outer;
        nonlinear_residual = residual.norm();
        if (iteration_log_)
            std::cout << "ECCAPFIM_ANDERSON step=" << current_timestep_ << " iteration=" << outer
                      << " residual=" << std::setprecision(17) << nonlinear_residual << std::endl;
        if (!residual.allFinite() || !std::isfinite(nonlinear_residual))
            throw std::runtime_error("ECCAPFIM: nonfinite coupled residual");
        if (nonlinear_residual <= coupled_atol_ + coupled_rtol_ * initial_norm) {
            converged = true;
            break;
        }
        if (outer == coupled_max_outer_)
            break;
        if (mass_matrix_every_iteration_)
            mass_matrix_ready_ = false;
        const Vector g = -apply_right_preconditioner(residual);
        Vector next = state + g;
        const int m = static_cast<int>(history_x.size());
        if (m > 0) {
            Eigen::MatrixXd dX(state.size(), m), dG(state.size(), m);
            for (int j = 0; j < m; ++j) {
                dX.col(j) = state - history_x[j];
                dG.col(j) = g - history_g[j];
            }
            const Vector gamma = dG.colPivHouseholderQr().solve(g);
            const Vector accelerated = next - (dX + dG) * gamma;
            if (accelerated.allFinite())
                next = accelerated;
        }
        history_x.push_back(state);
        history_g.push_back(g);
        if (static_cast<int>(history_x.size()) > coupled_history_) {
            history_x.erase(history_x.begin());
            history_g.erase(history_g.begin());
        }
        state = next;
        residual = coupled_residual(state);
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: coupled Anderson field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
}

void SimulationEccapfim::solve_fields_coupled() {
    RECORD_TIMER;
    Vector state = pack6(fieldE, fieldB);
    Vector residual = coupled_residual(state);
    if (coupled_solver_ == CoupledSolver::Anderson)
        solve_fields_anderson(state, residual);
    else if (coupled_solver_ == CoupledSolver::Ngmres)
        solve_fields_ngmres(state, residual);
    else
        solve_fields_newton_krylov(state, residual);
    // B_h is the midpoint field; promote both fields to integer time levels.
    unpack6(state, fieldEp, fieldBHalf);
    fieldEn = 2.0 * fieldEp - fieldE;
    fieldBn = 2.0 * fieldBHalf - fieldB;
    algorithmsECCAPFIM::copy_periodic_nodes(fieldEn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    algorithmsECCAPFIM::copy_periodic_nodes(fieldBn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
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

void SimulationEccapfim::make_diagnostic(int timestep) {
    RECORD_TIMER;
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
    write_diagnostics(timestep, diag,
                      {domain, species, fieldEn, fieldBn, fieldE_external, fieldBInit, bc_handler,
                       diagnostic_density_shape()},
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
