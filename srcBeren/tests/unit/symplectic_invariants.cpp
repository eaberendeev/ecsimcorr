#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <numbers>
#include <random>

#include "coupling/gather.h"
#include "diagnostics/DensityMeasurements.h"
#include "simulation/symplectic/simulation_symplectic.h"
#include "simulation/symplectic/whitney_pic.h"
#include "symplectic_reference.h"
#include "test_support.h"

namespace {
constexpr int g = GHOST_CELLS;

Vector3R sample_whitney_electric(const Domain& domain, const Field3d& field, const Vector3R& position) {
    const Vector3R coord = domain.to_cell_coordinates(position);
    const Stencil3D<2> stencil = make_linear_stencil_at(coord, checked_floor_index(coord));
    return gather_whitney_electric(field, stencil);
}

nlohmann::json settings(const Vector3I& cells, const Vector3R& h, double dt) {
    return {{"Scheme", "symplectic"},
            {"Dt", dt},
            {"Dx", h[0]},
            {"Dy", h[1]},
            {"Dz", h[2]},
            {"NumCellsX", cells[0]},
            {"NumCellsY", cells[1]},
            {"NumCellsZ", cells[2]},
            {"StartFromTime", 0},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

std::unique_ptr<SimulationSymplectic> make_simulation(const nlohmann::json& config, const nlohmann::json& species) {
    char name[] = "symplectic";
    char* argv[] = {name};
    std::unique_ptr<SimulationSymplectic> sim =
        std::make_unique<SimulationSymplectic>(config, nlohmann::json{{"particles", species}}, 1, argv);
    sim->init();
    return sim;
}

// Electrons with q/m = -1 and macro-weight `weight`; the box starts empty.
std::unique_ptr<SimulationSymplectic> make_beam(const Vector3I& cells, const Vector3R& h, double weight,
                                                double dt = 0.1) {
    return make_simulation(settings(cells, h, dt), nlohmann::json::array({check::species("Beam", -1, 1, weight, 1)}));
}

void fill_unique_nodes(const Domain& domain, Field3d& field, double amplitude, double frequency) {
    for (int i = 0; i < field.capacity(); ++i) field[i] = amplitude * std::sin(frequency * i);
    symplectic_reference::fill_periodic_ghosts(domain, field);
}

void clear(ParticlesArray& particles) {
    for (int cell = 0; cell < particles.size(); ++cell) particles(cell).clear();
}

std::vector<Particle> all_particles(const ParticlesArray& particles) {
    std::vector<Particle> result;
    for (int cell = 0; cell < particles.size(); ++cell)
        result.insert(result.end(), particles(cell).begin(), particles(cell).end());
    return result;
}

double max_abs(const Field3d& field) {
    return field.data().lpNorm<Eigen::Infinity>();
}

Field3d divergence(const SimulationSymplectic& sim, const Field3d& edges) {
    Field3d result(sim.domain.size(), 1);
    spmv(sim.maxwell.divE, edges, result);
    return result;
}

Field3d charge_density(const SimulationSymplectic& sim) {
    return measure_charge_density(sim.species, sim.domain, sim.bc_handler, ShapeType::Linear);
}

void check_periodic_position(const Domain& domain, const Vector3R& actual, const Vector3R& expected, double tolerance,
                             const std::string& message) {
    for (int d = 0; d < 3; ++d) {
        const double extent = domain.cell_size(d) * domain.num_cells(d);
        const double difference = actual[d] - expected[d];
        check::near(difference - extent * std::round(difference / extent), 0, tolerance, message);
    }
}

// Advance the bucketed particles and their independent copies, then compare the
// trajectories, the deposit and the local charge continuity of the production flow.
void compare_with_reference(SimulationSymplectic& sim, bool reverse, const std::string& label) {
    const Domain& domain = sim.domain;
    ParticlesArray& beam = *sim.species.at("Beam");
    const double qm = beam.charge() / beam.mass(), qw = beam.charge() * beam.mpw();
    const Field3d density_before = charge_density(sim);

    std::vector<Particle> expected = all_particles(beam);
    Field3d reference = Field3d::Zero(domain.size(), 3), unused = reference;
    for (Particle& particle : expected)
        for (int step = 0; step < 3; ++step)
            symplectic_reference::kinetic_flow_reference(domain, particle, reverse ? 2 - step : step, 1., qm, qw,
                                                         unused, sim.magnetic, reference);
    symplectic_reference::fill_periodic_ghosts(domain, reference);

    Field3d current = Field3d::Zero(domain.size(), 3);
    check::require(symplectic::kinetic_flow(beam, sim.magnetic, 1., reverse, current), label + ": admissible paths");
    sim.bc_handler.apply_to_fields(current, FieldType::CURRENT, domain);
    const std::vector<Particle> actual = all_particles(beam);
    for (size_t i = 0; i < actual.size(); ++i) {
        check_periodic_position(domain, actual[i].coord, expected[i].coord, 2e-13, label + ": trajectory");
        check::vector(actual[i].velocity, expected[i].velocity, 2e-13, label + ": magnetic path integral");
    }
    check::near((current - reference).norm(), 0, 3e-13, label + ": path deposit");

    sim.bc_handler.apply_to_particles(beam, sim.species, domain);
    beam.update_cells(domain);
    check::require(all_particles(beam).size() == actual.size(), label + ": periodic wrap keeps every particle");
    check::near(max_abs(charge_density(sim) - density_before + divergence(sim, current)), 0, 3e-14,
                label + ": local continuity across cells and periodic seams");
}

void add_axis_paths(ParticlesArray& beam, const Vector3R& h, int axis, const std::vector<double>& starts,
                    const std::vector<double>& travels) {
    for (double start : starts)
        for (double travel : travels) {
            Particle particle(Vector3R(0.31 * h[0], 0.57 * h[1], 0.71 * h[2]), Vector3R(0.13, -0.07, 0.09));
            particle.coord[axis] = start;
            particle.velocity[axis] = travel * h[axis];
            beam.add_particle(particle);
        }
}

void kinetic_flows() {
    for (const Vector3I& cells : {Vector3I(4, 3, 5), Vector3I(1, 2, 3), Vector3I(8, 7, 6)}) {
        const Vector3R h(0.7, 1.1, 0.9);
        std::unique_ptr<SimulationSymplectic> sim = make_beam(cells, h, 0.3);
        const Domain& domain = sim->domain;
        ParticlesArray& beam = *sim->species.at("Beam");
        const Vector3R extent(cells[0] * h[0], cells[1] * h[1], cells[2] * h[2]);

        for (bool reverse : {false, true}) {
            // Without B the velocities are constant: paths of up to exactly one cell
            // start on grid planes, at the seam and at the adjacent machine numbers.
            clear(beam);
            sim->magnetic.setZero();
            for (int axis = 0; axis < 3; ++axis) {
                add_axis_paths(beam, h, axis, {0., h[axis], std::nextafter(h[axis], 2 * h[axis]), 0.37 * h[axis]},
                               {-1., -0.999, 0., 0.999, 1.});
                // The last number below a plane may be binned on either side of it,
                // which leaves no room for a path of exactly one cell.
                add_axis_paths(beam, h, axis, {std::nextafter(h[axis], 0.), std::nextafter(extent[axis], 0.)},
                               {-0.999, 0., 0.999});
            }
            compare_with_reference(*sim, reverse, "one-cell paths");

            // B rotates the velocity between the subflows of one particle.
            clear(beam);
            fill_unique_nodes(domain, sim->magnetic, 0.05, 0.73);
            std::mt19937_64 rng(1729);
            std::uniform_real_distribution<double> unit(0, 1);
            for (int i = 0; i < 400; ++i) {
                Vector3R x, v;
                for (int d = 0; d < 3; ++d) {
                    x[d] = extent[d] * unit(rng);
                    v[d] = 1.2 * h[d] * (unit(rng) - 0.5);
                }
                beam.add_particle(Particle(x, v));
            }
            compare_with_reference(*sim, reverse, "magnetized paths");
        }

        Field3d electric = Field3d::Zero(domain.size(), 3);
        fill_unique_nodes(domain, electric, 1, 0.73);
        for (const Particle& particle : all_particles(beam))
            check::vector(sample_whitney_electric(domain, electric, particle.coord),
                          symplectic_reference::gather_reference(domain, electric, particle.coord, false), 2e-15,
                          "edge gather against tensor-product definition");
    }
}

void displacement_limit() {
    std::unique_ptr<SimulationSymplectic> sim = make_beam({3, 3, 3}, {0.7, 1.1, 0.9}, 0.3);
    ParticlesArray& beam = *sim->species.at("Beam");
    for (int axis = 0; axis < 3; ++axis)
        for (double travel : {-2.1, 2.1, std::numeric_limits<double>::quiet_NaN()}) {
            clear(beam);
            Particle particle({0.31, 0.41, 0.29}, {0, 0, 0});
            particle.velocity[axis] = travel * sim->domain.cell_size(axis);
            beam.add_particle(particle);
            Field3d current = Field3d::Zero(sim->domain.size(), 3);
            check::require(!symplectic::kinetic_flow(beam, sim->magnetic, 1., false, current),
                           "a path beyond the neighbor cell must be rejected");
            check::vector(all_particles(beam).front().coord, particle.coord, 0, "rejected path must not move");
            check::near(current.norm(), 0, 0, "rejected path must not deposit current");
        }

    bool rejected = false;
    try {
        sim->make_step(1);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    check::require(rejected, "the step must stop on an inadmissible path");
}

// The scheme state with a single macro-particle.
struct State {
    Field3d e, b;
    Particle p;
};

struct Model {
    std::unique_ptr<SimulationSymplectic> sim;

    Model(const Vector3I& cells, const Vector3R& h, double weight, double dt) : sim(make_beam(cells, h, weight, dt)) {
        load({sim->electric, sim->magnetic, Particle(Vector3R(0.37, 0.41, 0.29), Vector3R(0.13, -0.07, 0.09))});
    }

    ParticlesArray& beam() const {
        return *sim->species.at("Beam");
    }

    void load(const State& state) {
        sim->electric = state.e;
        sim->magnetic = state.b;
        clear(beam());
        beam().add_particle(state.p);
    }

    State save() const {
        return {sim->electric, sim->magnetic, all_particles(beam()).front()};
    }

    void step(double dt) {
        sim->dt_ = dt;
        sim->prepare_step(1);
        sim->make_step(1);
    }

    double energy() const {
        const Vector3R v = all_particles(beam()).front().velocity;
        return sim->field_energy(sim->electric) + sim->field_energy(sim->magnetic) +
               get_energy_particle(v, beam().mass(), beam().mpw());
    }
};

void exact_subflows() {
    Model model({4, 3, 5}, {0.7, 1.1, 0.9}, 0.3, 0.07);
    SimulationSymplectic& sim = *model.sim;
    const Domain& domain = sim.domain;

    // In the reversed order only H_x moves a particle that starts with v = (vx, 0, 0).
    State state = model.save();
    state.p.velocity = {0.13, 0, 0};
    blas::fill(state.b, 0.6);
    model.load(state);
    Field3d current = Field3d::Zero(domain.size(), 3);
    check::require(symplectic::kinetic_flow(model.beam(), sim.magnetic, 0.4, true, current), "admissible path");
    check::vector(model.save().p.velocity,
                  state.p.velocity + Vector3R(0.4 * 0.13, 0, 0).cross(Vector3R(0.6, 0.6, 0.6)) * -1.0, 1e-15,
                  "exact magnetic line integral and Lorentz sign");

    // Discrete curl and Whitney reconstruction commute inside a cell.
    fill_unique_nodes(domain, state.e, 1, 0.73);
    const Field3d curl = sim.maxwell.curlE * state.e;
    const Vector3R x = state.p.coord;
    Vector3R finite_curl(0, 0, 0);
    for (int d = 0; d < 3; ++d) {
        const int a = (d + 1) % 3, b = (d + 2) % 3;
        Vector3R plus = x, minus = x;
        plus[a] += 1e-6;
        minus[a] -= 1e-6;
        finite_curl[d] =
            (sample_whitney_electric(domain, state.e, plus)[b] - sample_whitney_electric(domain, state.e, minus)[b]) /
            2e-6;
        plus = x;
        minus = x;
        plus[b] += 1e-6;
        minus[b] -= 1e-6;
        finite_curl[d] -=
            (sample_whitney_electric(domain, state.e, plus)[a] - sample_whitney_electric(domain, state.e, minus)[a]) /
            2e-6;
    }
    check::vector(symplectic_reference::gather_reference(domain, curl, x, true), finite_curl, 5e-10,
                  "Whitney curl compatibility");
    check::near(symplectic::max_magnetic_divergence(domain, curl), 0, 3e-15, "div curl E = 0");

    state.b = curl;
    state.p.velocity = {0.13, -0.07, 0.09};
    model.load(state);
    model.step(0.07);
    model.step(-0.07);
    const State reversed = model.save();
    check::vector(reversed.p.coord, state.p.coord, 1e-14, "time reversal position");
    check::vector(reversed.p.velocity, state.p.velocity, 1e-14, "time reversal velocity");
    check::near((reversed.e - state.e).norm(), 0, 3e-14, "time reversal E");
    check::near((reversed.b - state.b).norm(), 0, 3e-14, "time reversal B");
}

// Uniform cold electron plasma with immobile uniform positive background:
// dE/dt = v, dv/dt = -E, omega_p=1. One periodic cell gives uniform Whitney fields.
double plasma_error(double dt, int steps, bool long_run = false) {
    Model model({1, 1, 1}, {1, 1, 1}, 1, dt);
    State state = model.save();
    state.p.velocity = {0.1, 0, 0};
    model.load(state);
    double maximum = 0;
    for (int i = 0; i < steps; ++i) {
        model.step(dt);
        const double energy = model.energy();
        check::require(std::isfinite(energy), "non-finite cold plasma energy");
        maximum = std::max(maximum, std::abs(energy / 0.005 - 1));
    }
    if (long_run) {
        check::require(maximum < 0.0007, "cold plasma long-time energy envelope");
        std::cout << "cold plasma: " << steps << " steps, max relative energy error=" << maximum << '\n';
    }
    const double t = steps * dt;
    const State last = model.save();
    return std::hypot(last.e(g, g, g, 0) - 0.1 * std::sin(t), last.p.velocity[0] - 0.1 * std::cos(t));
}

// The exactly soluble split stages do not make their composition unconditionally
// stable. In a homogeneous electron plasma the map has trace 2-(omega_p*dt)^2.
void plasma_stability_limit() {
    for (double dt : {1.5, 1.9, 2.1}) {
        Model model({1, 1, 1}, {10, 10, 10}, 1, dt);
        State state = model.save();
        state.p.velocity = {1e-8, 0, 0};
        model.load(state);
        double maximum = 0;
        for (int i = 0; i < 25; ++i) {
            model.step(dt);
            const double ratio = model.energy() / 5e-17;
            check::require(std::isfinite(ratio), "non-finite plasma stability experiment");
            maximum = std::max(maximum, ratio);
        }
        if (dt < 2)
            check::require(maximum < 1.01, "resolved stable plasma map at 1 < omega_p*dt < 2");
        else
            check::require(maximum > 1e10, "plasma map must expose instability beyond omega_p*dt=2");
        std::cout << "omega_p*dt=" << dt << ": maximum H/H0=" << maximum << '\n';
    }
}

double vacuum_error(double dt, int steps) {
    std::unique_ptr<SimulationSymplectic> sim = make_beam({8, 1, 1}, {0.7, 1.1, 0.9}, 1, dt);
    for (int i = 0; i < 8; ++i) sim->electric(i + g, g, g, 1) = std::cos(2 * std::numbers::pi * i / 8);
    symplectic_reference::fill_periodic_ghosts(sim->domain, sim->electric);
    const Field3d initial = sim->electric;
    const Field3d curl = sim->maxwell.curlE * initial;
    for (int i = 1; i <= steps; ++i) sim->make_step(i);
    const double omega = 2 * std::sin(std::numbers::pi / 8) / 0.7;
    const double t = steps * dt;
    return (sim->electric - std::cos(omega * t) * initial).norm() +
           (sim->magnetic + std::sin(omega * t) / omega * curl).norm();
}

// Independent noncanonical Poisson-map check, including particle-field blocks.
// z=(E,B,x,v), H=(sum E^2+B^2 + w*m*v^2)/2, m=1, q=-1, w=0.2.
// A symplectic splitting must satisfy DF J(z) DF^T = J(F(z)) on div B=0.
void poisson_map() {
    constexpr int nf = 24, nx = 48, nv = 51, total = 54;
    constexpr double dt = 0.04;
    Model model({2, 2, 2}, {0.7, 0.9, 1.1}, 0.2, dt);
    const SimulationSymplectic& sim = *model.sim;
    const Domain& domain = sim.domain;

    State initial = model.save();
    fill_unique_nodes(domain, initial.e, 0.03, 0.71);
    initial.b = -1.0 * (sim.maxwell.curlE * initial.e);
    for (int i = 0; i < initial.b.capacity(); ++i) initial.b[i] += 0.13;

    std::function<Eigen::VectorXd(const State&)> pack = [&](const State& state) {
        Eigen::VectorXd z(total);
        int r = 0;
        for (int i = g; i < g + 2; ++i)
            for (int j = g; j < g + 2; ++j)
                for (int k = g; k < g + 2; ++k)
                    for (int d = 0; d < 3; ++d, ++r) {
                        z[r] = state.e(i, j, k, d);
                        z[nf + r] = state.b(i, j, k, d);
                    }
        for (int d = 0; d < 3; ++d) {
            z[nx + d] = state.p.coord[d];
            z[nv + d] = state.p.velocity[d];
        }
        return z;
    };
    std::function<State(const Eigen::VectorXd&)> unpack = [&](const Eigen::VectorXd& z) {
        State state{Field3d::Zero(domain.size(), 3), Field3d::Zero(domain.size(), 3), initial.p};
        int r = 0;
        for (int i = g; i < g + 2; ++i)
            for (int j = g; j < g + 2; ++j)
                for (int k = g; k < g + 2; ++k)
                    for (int d = 0; d < 3; ++d, ++r) {
                        state.e(i, j, k, d) = z[r];
                        state.b(i, j, k, d) = z[nf + r];
                    }
        for (int d = 0; d < 3; ++d) {
            state.p.coord[d] = z[nx + d];
            state.p.velocity[d] = z[nv + d];
        }
        symplectic_reference::fill_periodic_ghosts(domain, state.e);
        symplectic_reference::fill_periodic_ghosts(domain, state.b);
        return state;
    };
    std::function<State(const State&)> advance = [&](const State& state) {
        model.load(state);
        model.step(dt);
        return model.save();
    };
    std::function<Eigen::MatrixXd(const State&)> bracket = [&](const State& state) {
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(total, total);
        for (int r = 0; r < nf; ++r) {
            Eigen::VectorXd basis = Eigen::VectorXd::Zero(total);
            basis[r] = 1;
            State unit = unpack(basis);
            unit.b = -1.0 * (sim.maxwell.curlE * unit.e);
            const Eigen::VectorXd column = pack(unit);
            for (int t = 0; t < nf; ++t) {
                J(nf + t, r) = column[nf + t];
                J(r, nf + t) = -column[nf + t];
            }
            const Vector3R w = sample_whitney_electric(domain, unit.e, state.p.coord);
            for (int d = 0; d < 3; ++d) {
                J(nv + d, r) = -w[d];
                J(r, nv + d) = w[d];
            }
        }
        for (int d = 0; d < 3; ++d) {
            J(nx + d, nv + d) = 5;
            J(nv + d, nx + d) = -5;
        }
        const Vector3R b = symplectic_reference::gather_reference(domain, state.b, state.p.coord, true);
        for (int d = 0; d < 3; ++d) {
            const int a = (d + 1) % 3, c = (d + 2) % 3;
            J(nv + a, nv + c) = -5 * b[d];
            J(nv + c, nv + a) = 5 * b[d];
        }
        return J;
    };

    const Eigen::VectorXd z = pack(initial);
    Eigen::MatrixXd derivative(total, total);
    constexpr double epsilon = 1e-6;
    for (int r = 0; r < total; ++r) {
        Eigen::VectorXd plus = z, minus = z;
        plus[r] += epsilon;
        minus[r] -= epsilon;
        derivative.col(r) = (pack(advance(unpack(plus))) - pack(advance(unpack(minus)))) / (2 * epsilon);
    }
    const Eigen::MatrixXd defect = derivative * bracket(initial) * derivative.transpose() - bracket(advance(initial));
    check::near(defect.cwiseAbs().maxCoeff(), 0, 3e-8, "full particle-field noncanonical Poisson map");
    std::cout << "Poisson-map maximum defect=" << defect.cwiseAbs().maxCoeff() << '\n';
}

void convergence() {
    const double p1 = plasma_error(0.1, 10), p2 = plasma_error(0.05, 20), p3 = plasma_error(0.025, 40);
    const double v1 = vacuum_error(0.1, 10), v2 = vacuum_error(0.05, 20);
    check::require(p1 / p2 > 3.9 && p1 / p2 < 4.1 && p2 / p3 > 3.9 && p2 / p3 < 4.1,
                   "second-order convergence of cold plasma frequency and phase");
    check::require(v1 / v2 > 3.9 && v1 / v2 < 4.1, "second-order vacuum Maxwell wave");
    std::cout << "convergence ratios: plasma=" << p1 / p2 << ", vacuum=" << v1 / v2 << '\n';
    plasma_error(0.05, 20000, true);
}

bool lexicographic(const Particle& a, const Particle& b) {
    return std::lexicographical_compare(a.coord.data, a.coord.data + 3, b.coord.data, b.coord.data + 3);
}

void parallel_deposits() {
    const int original_threads = omp_get_max_threads();
    std::function<std::unique_ptr<SimulationSymplectic>(int)> run = [&](int threads) {
        // The Yee operators are partitioned for the thread count of their construction.
        omp_set_num_threads(threads);
        std::unique_ptr<SimulationSymplectic> sim = make_simulation(
            settings({24, 3, 3}, {1, 1, 1}, 0.1),
            nlohmann::json::array({check::species("Beam", -1, 1, 1, 512), check::species("Ions", 1, 25, 1, 512),
                                   check::species("Gas", 0, 2, 1, 512)}));
        std::mt19937_64 rng(71231);
        std::uniform_real_distribution<double> unit(0, 1);
        for (const std::string name : {"Beam", "Ions", "Gas"}) {
            ParticlesArray& sp = *sim->species.at(name);
            const int count = sp.is_neutral() ? 5000 : 10000;
            for (int i = 0; i < count; ++i) {
                Vector3R x(24 * unit(rng), 3 * unit(rng), 3 * unit(rng));
                Vector3R v(0.2 * (unit(rng) - 0.5), 0.2 * (unit(rng) - 0.5), 0.2 * (unit(rng) - 0.5));
                if (i % 31 == 0) {
                    x[0] = 23.999;
                    v[0] = 0.1;
                }
                // Neutral crosses the periodic seam without depositing current.
                if (sp.is_neutral() && i == 0) {
                    x = {23.99, 1.2, 0.4};
                    v = {1, 0, 0};
                }
                sp.add_particle(Particle(x, v));
            }
        }

        const Field3d initial_defect = sim->gauss_defect();
        Field3d old_density = charge_density(*sim);
        Field3d current(sim->domain.size(), 3);
        for (int step = 1; step <= 3; ++step) {
            sim->prepare_step(step);
            sim->make_step(step);
            check::near(sim->species.at("Gas")->currentOnGrid.norm(), 0, 0, "neutral species must not deposit current");
            const Field3d density = charge_density(*sim);
            sim->collect_current(current);
            check::near(max_abs((density - old_density) * (1 / 0.1) + divergence(*sim, current)), 0, 3e-13,
                        "parallel local charge continuity");
            check::near(max_abs(sim->gauss_defect() - initial_defect), 0, 3e-13, "parallel Gauss preservation");
            old_density = density;
        }
        return sim;
    };

    const std::unique_ptr<SimulationSymplectic> serial = run(1), parallel = run(4);
    check::near((serial->electric - parallel->electric).norm(), 0, 2e-13, "thread-count electric agreement");
    check::near((serial->magnetic - parallel->magnetic).norm(), 0, 2e-13, "thread-count magnetic agreement");
    for (const Species::value_type& entry : serial->species) {
        const ParticlesArray& sa = *entry.second;
        const ParticlesArray& sb = *parallel->species.at(entry.first);
        check::near((sa.currentOnGrid - sb.currentOnGrid).norm(), 0, 2e-13, "thread-count current agreement");
        for (int cell = 0; cell < sa.size(); ++cell) {
            // The order inside a bucket depends on the thread schedule.
            std::vector<Particle> a = sa(cell), b = sb(cell);
            check::require(a.size() == b.size(), "thread-count bucket population");
            std::sort(a.begin(), a.end(), lexicographic);
            std::sort(b.begin(), b.end(), lexicographic);
            for (size_t i = 0; i < a.size(); ++i) {
                check::vector(a[i].coord, b[i].coord, 3e-14, "thread-count position");
                check::vector(a[i].velocity, b[i].velocity, 3e-14, "thread-count velocity");
                if (sa.is_neutral() && a[i].velocity[0] == 1)
                    check::vector(a[i].coord, {0.29, 1.2, 0.4}, 2e-14, "neutral periodic rebin");
            }
        }
    }
    omp_set_num_threads(original_threads);
}

void simulation() {
    nlohmann::json cfg = settings({3, 3, 3}, {1, 1, 1}, 0.1);
    cfg["diagnostics"]["outputs"] = {{{"type", "energy_balance"}, {"interval", 1}},
                                     {{"type", "recovery"}, {"interval", 1}}};
    const nlohmann::json particles =
        nlohmann::json::array({check::species("Beam", -1, 1, 1, 1), check::species("Gas", 0, 2, 1, 1)});
    std::unique_ptr<SimulationSymplectic> owner = make_simulation(cfg, particles);
    SimulationSymplectic& sim = *owner;
    for (const Species::value_type& entry : sim.species)
        check::near(entry.second->currentOnGrid.norm(), 0, 0, "step-zero current must be initialized");

    // A coordinate just below the box rounds to its excluded upper bound.
    check::require(sim.domain.contains(sim.bc_handler.wrap_periodic({-1e-17, 1.5, 3 - 1e-17}, sim.domain)),
                   "periodic wrap must return a position inside the box");

    ParticlesArray& beam = *sim.species.at("Beam");
    beam.add_particle(Particle({2.99, 0.8, 1.2}, {0.4, -0.1, 0.03}));
    ParticlesArray& gas = *sim.species.at("Gas");
    gas.add_particle(Particle({0.02, 1.2, 0.6}, {-1, 0.03, 0.02}));
    const Field3d rho_old = charge_density(sim);
    const Field3d defect_old = sim.gauss_defect();
    sim.prepare_step(1);
    sim.make_step(1);
    check::near(max_abs((charge_density(sim) - rho_old) * (1 / 0.1) + divergence(sim, beam.currentOnGrid)), 0, 2e-14,
                "complete Simulation step charge continuity");
    check::near(max_abs(sim.gauss_defect() - defect_old), 0, 3e-15, "Gauss defect preserved");
    check::near(symplectic::max_magnetic_divergence(sim.domain, sim.magnetic), 0, 1e-15, "div B preserved");
    for (const Particle& p : all_particles(gas)) {
        check::vector(p.coord, {2.92, 1.203, 0.602}, 1e-14, "neutral periodic drift");
        check::vector(p.velocity, {-1, 0.03, 0.02}, 0, "neutral velocity");
    }
    for (const Species::value_type& entry : sim.species) {
        const ParticlesArray& sp = *entry.second;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const Particle& p : sp(cell)) {
                const Vector3I index = sim.domain.storage_cell_index(p.coord);
                check::require(sim.domain.sind(index[0], index[1], index[2]) == cell,
                               "particle not rebinned into its final cell");
            }
    }

    sim.make_diagnostic(0);
    sim.make_diagnostic(1);
    sim.finalize();
    sim.finalize();
    std::ifstream table("energy.txt");
    std::string header, row;
    std::getline(table, header);
    check::require(header.find("gaussDrift") != std::string::npos && header.find("energyError") != std::string::npos,
                   "symplectic energy output missing structural diagnostics");
    int rows = 0;
    while (std::getline(table, row))
        if (!row.empty())
            ++rows;
    check::require(rows == 2, "energy output step zero and step one");

    // Recovery stores synchronous E,B,x,v. Compare a continued step with an
    // actual disk restart, including a charged particle that crossed the seam.
    sim.prepare_step(2);
    sim.make_step(2);
    const std::filesystem::path original_path = std::filesystem::current_path();
    std::filesystem::create_directory("restart");
    std::filesystem::current_path("restart");
    nlohmann::json restart_cfg = cfg;
    restart_cfg["StartFromTime"] = 1;
    restart_cfg["StartTimeStep"] = 1;
    restart_cfg["diagnostics"]["outputs"] = nlohmann::json::array();
    std::unique_ptr<SimulationSymplectic> restart = make_simulation(restart_cfg, particles);
    restart->prepare_step(2);
    restart->make_step(2);
    check::near((restart->electric - sim.electric).norm(), 0, 1e-14, "restart electric field");
    check::near((restart->magnetic - sim.magnetic).norm(), 0, 1e-14, "restart magnetic field");
    for (const Species::value_type& entry : sim.species) {
        const std::vector<Particle> continued = all_particles(*entry.second);
        const std::vector<Particle> restored = all_particles(*restart->species.at(entry.first));
        check::require(continued.size() == restored.size(), "restart population");
        for (size_t p = 0; p < continued.size(); ++p) {
            check::vector(restored[p].coord, continued[p].coord, 1e-14, "restart position");
            check::vector(restored[p].velocity, continued[p].velocity, 1e-14, "restart velocity");
        }
    }
    restart->make_diagnostic(0);
    restart->finalize();
    std::filesystem::current_path(original_path);

    for (int variant = 0; variant < 8; ++variant) {
        nlohmann::json bad = cfg;
        if (variant == 0)
            bad["Dt"] = 1;
        if (variant == 1)
            bad["Boundary_conditions"] = nlohmann::json::array();
        if (variant == 2)
            bad["Collider"] = "Coulomb";
        if (variant == 3)
            bad["DampingType"] = "CircleXY";
        if (variant == 4)
            bad["ExternalFieldB"] = {1, 2, 3};
        if (variant == 5)
            bad["Symplectic"] = {{"order", 4}};
        if (variant == 6)
            bad["Dt"] = std::numeric_limits<double>::quiet_NaN();
        if (variant == 7)
            bad["Symplectic"] = {{"initial_magnetic", {0, 0, std::numeric_limits<double>::infinity()}}};
        bool rejected = false;
        try {
            make_simulation(bad, particles);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check::require(rejected, "unsupported symplectic configuration accepted");
    }
}
}   // namespace

int main() {
    return check::run([] {
        kinetic_flows();
        displacement_limit();
        exact_subflows();
        poisson_map();
        convergence();
        plasma_stability_limit();
        simulation();
        parallel_deposits();
    });
}
