#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <numbers>
#include <random>

#include "simulation/symplectic/simulation_symplectic.h"
#include "symplectic_reference.h"
#include "test_support.h"

namespace {
void kinetic_flow(const Domain& domain, Particle& particle, int axis, double dt, double qm, double qw,
                  Field3d& electric, const Field3d& magnetic, Field3d& integral) {
    Field3d delta = Field3d::Zero(domain.size(), 3);
    symplectic::advance_axis(domain, particle, axis, dt, qm, qw, magnetic, delta);
    electric.data() -= delta.data();
    integral.data() += delta.data();
}

struct State {
    Field3d e, b;
    Particle p;
};

State state(const Domain& domain) {
    return {Field3d::Zero(domain.size(), 3), Field3d::Zero(domain.size(), 3),
            Particle(Vector3R(0.37, 0.41, 0.29), Vector3R(0.13, -0.07, 0.09))};
}

void step(const Domain& domain, State& s, double dt, double qm = -1, double qw = -0.2) {
    const Domain& pic_domain = domain;
    Field3d integral = Field3d::Zero(domain.size(), 3);
    std::function<void(int, double)> drift = [&](int axis, double h) {
        kinetic_flow(pic_domain, s.p, axis, h, qm, qw, s.e, s.b, integral);
    };
    drift(0, dt / 2);
    drift(1, dt / 2);
    drift(2, dt / 2);
    symplectic::magnetic_flow(pic_domain, s.e, s.b, dt / 2);
    symplectic::electric_flow(pic_domain, s.b, s.e, dt);
    s.p.velocity += dt * qm * symplectic::gather<false>(pic_domain, s.e, s.p.coord);
    symplectic::magnetic_flow(pic_domain, s.e, s.b, dt / 2);
    drift(2, dt / 2);
    drift(1, dt / 2);
    drift(0, dt / 2);
    symplectic::fill_periodic_ghosts(pic_domain, s.e);
    symplectic::fill_periodic_ghosts(pic_domain, s.b);
}

void optimized_kernels() {
    Domain domain;
    domain.init({4, 3, 5}, {0.7, 1.1, 0.9});
    const Domain& pic_domain = domain;
    Field3d field = Field3d::Zero(domain.size(), 3);
    for (int i = 0; i < field.capacity(); ++i) field[i] = std::sin(0.73 * i);
    symplectic::fill_periodic_ghosts(pic_domain, field);
    const std::vector<Vector3R> positions = {{0, 0, 0},         {0.7, 1.1, 0.9},   {2.79999, 3.29999, 4.49999},
                                             {-0.01, 0.1, 4.6}, {7.2, -5.1, 13.4}, {0.31, 1.57, 2.71}};
    for (const Vector3R& x : positions) {
        for (bool magnetic : {false, true})
            check::vector((magnetic ? symplectic::gather<true>(pic_domain, field, x)
                                    : symplectic::gather<false>(pic_domain, field, x)),
                          symplectic_reference::gather_reference(pic_domain, field, x, magnetic), 2e-15,
                          "optimized Whitney gather against tensor-product definition");
        Field3d rho = Field3d::Zero(domain.size(), 1), reference = rho;
        symplectic::deposit_charge(pic_domain, rho, x, -0.3);
        symplectic_reference::deposit_charge_reference(pic_domain, reference, x, -0.3);
        check::near((rho.data() - reference.data()).norm(), 0, 2e-15, "optimized nodal charge");
        for (int axis = 0; axis < 3; ++axis)
            for (double velocity : {-19.3, 0., 21.2}) {
                Particle p(x, {0.13, -0.07, 0.09});
                p.velocity[axis] = velocity;
                Particle original = p;
                Field3d e = field, e_reference = field;
                Field3d integral = Field3d::Zero(domain.size(), 3), j_reference = integral;
                kinetic_flow(pic_domain, p, axis, 1.3, -1, -0.3, e, field, integral);
                symplectic_reference::kinetic_flow_reference(pic_domain, original, axis, 1.3, -1, -0.3, e_reference,
                                                             field, j_reference);
                check::vector(p.coord, original.coord, 2e-14, "optimized axis path endpoint");
                check::vector(p.velocity, original.velocity, 2e-13, "optimized magnetic path integral");
                check::near((e.data() - e_reference.data()).norm(), 0, 3e-14, "optimized electric response");
                check::near((integral.data() - j_reference.data()).norm(), 0, 3e-14, "optimized path deposit");
            }
    }
}

void bucket_flows() {
    for (const Vector3I& cells : {Vector3I(1, 2, 3), Vector3I(8, 7, 6)}) {
        Domain domain;
        domain.init(cells, {0.7, 1.1, 0.9});
        Field3d magnetic = Field3d::Zero(domain.size(), 3);
        for (int i = 0; i < magnetic.capacity(); ++i) magnetic[i] = 0.17 * std::sin(0.73 * i);
        symplectic::fill_periodic_ghosts(domain, magnetic);
        const Vector3R origin(0.7 * (cells[0] - 0.1), 0.03, 0.41);
        const Vector3I cell = domain.storage_cell_index(origin);
        const int bucket = domain.sind(cell[0], cell[1], cell[2]);
        for (bool reverse : {false, true}) {
            for (double qm : {0., -1.}) {
                std::vector<Particle> particles;
                for (int i = 0; i < 32; ++i)
                    particles.emplace_back(origin, Vector3R(0.13, -0.07, 0.09) * static_cast<double>(i - 16));
                // Exercise paths beyond the cached patch and through multiple periodic seams.
                particles.emplace_back(origin, Vector3R(-19.3, 21.2, -9.1));
                std::vector<Particle> expected = particles;
                Field3d integral = Field3d::Zero(domain.size(), 3), reference = integral, electric = integral;
                const double qw = qm == 0 ? 0 : -0.3;
                for (Particle& particle : expected)
                    for (int step = 0; step < 3; ++step) {
                        const int axis = reverse ? 2 - step : step;
                        symplectic_reference::kinetic_flow_reference(domain, particle, axis, 0.7, qm, qw, electric,
                                                                     magnetic, reference);
                    }
                symplectic::advance_bucket(domain, particles, bucket, reverse, 0.7, qm, qw, magnetic, integral);
                for (size_t i = 0; i < particles.size(); ++i) {
                    check::vector(particles[i].coord, expected[i].coord, 2e-13, "bucket trajectory against reference");
                    check::vector(particles[i].velocity, expected[i].velocity, 2e-13,
                                  "bucket velocity against reference");
                }
                check::near((integral.data() - reference.data()).norm(), 0, 3e-13, "bucket current against reference");
            }
        }
    }
}

void periodic_ghosts() {
    Domain domain;
    domain.init({1, 3, 2}, {0.7, 1.1, 0.9});
    Field3d field = Field3d::Zero(domain.size(), 3);
    const int g = GHOST_CELLS;
    for (int i = 0; i < 1; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 2; ++k)
                for (int d = 0; d < 3; ++d) field(i + g, j + g, k + g, d) = 100 * j + 10 * k + d;
    const int original_threads = omp_get_max_threads();
    omp_set_num_threads(1);
    const double serial_energy = symplectic::field_energy(domain, field);
    omp_set_num_threads(4);
    check::near(symplectic::field_energy(domain, field), serial_energy, 0, "field energy has fixed summation order");
    omp_set_num_threads(original_threads);
    const Field3d initial = field;
    symplectic::fill_periodic_ghosts(domain, field);
    const Vector3I size = domain.size();
    for (int i = 0; i < size[0]; ++i)
        for (int j = 0; j < size[1]; ++j)
            for (int k = 0; k < size[2]; ++k)
                for (int d = 0; d < 3; ++d)
                    check::near(field(i, j, k, d),
                                symplectic_reference::value(domain, initial, {i - g, j - g, k - g}, d), 0,
                                "periodic copies read only unique owned nodes");
}

void continuity() {
    Domain domain;
    domain.init({4, 3, 5}, {0.7, 1.1, 0.9});
    const Domain& pic_domain = domain;
    for (int axis = 0; axis < 3; ++axis)
        for (double velocity : {-19.3, -1.4, 0., 1.7, 21.2}) {
            State s = state(domain);
            // Start exactly on a grid plane to exercise the negative path convention.
            s.p.coord = {0.7, 1.1, 0.9};
            s.p.velocity[axis] = velocity;
            const Vector3R start = s.p.coord;
            Field3d old = Field3d::Zero(domain.size(), 1), rho = old;
            Field3d integral = Field3d::Zero(domain.size(), 3);
            symplectic::deposit_charge(pic_domain, old, start, -0.3);
            kinetic_flow(pic_domain, s.p, axis, 1.3, -1, -0.3, s.e, s.b, integral);
            symplectic::deposit_charge(pic_domain, rho, s.p.coord, -0.3);
            const int g = GHOST_CELLS;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 5; ++k) {
                        const double change = rho(i + g, j + g, k + g, 0) - old(i + g, j + g, k + g, 0);
                        check::near(change + symplectic::divergence(pic_domain, integral, i, j, k, false), 0, 3e-14,
                                    "local continuity across cells and periodic wraps");
                    }
        }
}

void exact_subflows() {
    Domain domain;
    domain.init({4, 3, 5}, {0.7, 1.1, 0.9});
    const Domain& pic_domain = domain;
    State s = state(domain);
    s.b.data().setConstant(0.6);
    const Particle start = s.p;
    Field3d integral = Field3d::Zero(domain.size(), 3);
    kinetic_flow(pic_domain, s.p, 0, 0.4, -2, 0, s.e, s.b, integral);
    check::vector(s.p.velocity,
                  start.velocity + Vector3R(0.4 * start.velocity[0], 0, 0).cross(Vector3R(0.6, 0.6, 0.6)) * -2.0, 1e-15,
                  "exact magnetic line integral and Lorentz sign");
    // Discrete curl and Whitney reconstruction commute inside a cell.
    for (int i = 0; i < s.e.capacity(); ++i) s.e[i] = std::sin(0.73 * i);
    symplectic::fill_periodic_ghosts(pic_domain, s.e);
    Field3d curl = Field3d::Zero(domain.size(), 3);
    symplectic::electric_flow(pic_domain, curl, s.e, -1);
    const Vector3R reconstructed = symplectic::gather<true>(pic_domain, curl, start.coord);
    Vector3R finite_curl(0, 0, 0);
    for (int d = 0; d < 3; ++d) {
        const int a = (d + 1) % 3, b = (d + 2) % 3;
        Vector3R plus = start.coord, minus = start.coord;
        plus[a] += 1e-6;
        minus[a] -= 1e-6;
        finite_curl[d] = (symplectic::gather<false>(pic_domain, s.e, plus)[b] -
                          symplectic::gather<false>(pic_domain, s.e, minus)[b]) /
                         2e-6;
        plus = start.coord;
        minus = start.coord;
        plus[b] += 1e-6;
        minus[b] -= 1e-6;
        finite_curl[d] -= (symplectic::gather<false>(pic_domain, s.e, plus)[a] -
                           symplectic::gather<false>(pic_domain, s.e, minus)[a]) /
                          2e-6;
    }
    check::vector(reconstructed, finite_curl, 5e-10, "Whitney curl compatibility");
    s.b = curl;
    symplectic::fill_periodic_ghosts(pic_domain, s.b);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 5; ++k)
                check::near(symplectic::divergence(pic_domain, s.b, i, j, k, true), 0, 3e-15, "div curl E = 0");
    const State original = s;
    step(domain, s, 0.07);
    step(domain, s, -0.07);
    check::vector(s.p.coord, original.p.coord, 1e-14, "time reversal position");
    check::vector(s.p.velocity, original.p.velocity, 1e-14, "time reversal velocity");
    check::near((s.e.data() - original.e.data()).norm(), 0, 3e-14, "time reversal E");
    check::near((s.b.data() - original.b.data()).norm(), 0, 3e-14, "time reversal B");
}

// Uniform cold electron plasma with immobile uniform positive background:
// dE/dt = v, dv/dt = -E, omega_p=1. One periodic cell gives uniform Whitney fields.
double plasma_error(double dt, int steps, bool long_run = false) {
    Domain domain;
    domain.init({1, 1, 1}, {1, 1, 1});
    const Domain& pic_domain = domain;
    State s = state(domain);
    s.p.velocity = {0.1, 0, 0};
    double maximum = 0;
    for (int i = 0; i < steps; ++i) {
        step(domain, s, dt, -1, -1);
        const double energy = symplectic::field_energy(pic_domain, s.e) + 0.5 * s.p.velocity.dot(s.p.velocity);
        check::require(std::isfinite(energy), "non-finite cold plasma energy");
        maximum = std::max(maximum, std::abs(energy / 0.005 - 1));
    }
    if (long_run) {
        check::require(maximum < 0.0007, "cold plasma long-time energy envelope");
        std::cout << "cold plasma: 20000 steps, max relative energy error=" << maximum << '\n';
    }
    const double t = steps * dt;
    return std::hypot(s.e(1, 1, 1, 0) - 0.1 * std::sin(t), s.p.velocity[0] - 0.1 * std::cos(t));
}

// The exactly soluble split stages do not make their composition unconditionally
// stable. In a homogeneous electron plasma the map has trace 2-(omega_p*dt)^2.
void plasma_stability_limit() {
    Domain domain;
    domain.init({1, 1, 1}, {10, 10, 10});
    const Domain& pic_domain = domain;
    for (double dt : {1.5, 1.9, 2.1}) {
        State s = state(domain);
        s.p.velocity = {1e-8, 0, 0};
        double maximum = 0;
        for (int i = 0; i < 25; ++i) {
            step(domain, s, dt, -1, -1);
            const double ratio =
                (symplectic::field_energy(pic_domain, s.e) + 0.5 * s.p.velocity.dot(s.p.velocity)) / 5e-17;
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
    Domain domain;
    domain.init({8, 1, 1}, {0.7, 1.1, 0.9});
    State s = state(domain);
    const Domain& pic_domain = domain;
    for (int i = 0; i < 8; ++i) s.e(i + 1, 1, 1, 1) = std::cos(2 * std::numbers::pi * i / 8);
    symplectic::fill_periodic_ghosts(pic_domain, s.e);
    Field3d initial = s.e;
    Field3d curl = Field3d::Zero(domain.size(), 3);
    symplectic::electric_flow(pic_domain, curl, initial, -1);
    symplectic::fill_periodic_ghosts(pic_domain, curl);
    for (int i = 0; i < steps; ++i) step(domain, s, dt, 0, 0);
    const double omega = 2 * std::sin(std::numbers::pi / 8) / 0.7;
    const double t = steps * dt;
    return (s.e.data() - std::cos(omega * t) * initial.data()).norm() +
           (s.b.data() + std::sin(omega * t) / omega * curl.data()).norm();
}

// Independent noncanonical Poisson-map check, including particle-field blocks.
// z=(E,B,x,v), H=(sum E^2+B^2 + w*m*v^2)/2, m=1, q=-1, w=0.2.
// A symplectic splitting must satisfy DF J(z) DF^T = J(F(z)) on div B=0.
void poisson_map() {
    Domain domain;
    domain.init({2, 2, 2}, {0.7, 0.9, 1.1});
    const Domain& pic_domain = domain;
    constexpr int nf = 24, nx = 48, nv = 51, total = 54;
    State s = state(domain);
    for (int i = 0; i < s.e.capacity(); ++i) s.e[i] = 0.03 * std::sin(0.71 * i);
    symplectic::fill_periodic_ghosts(pic_domain, s.e);
    symplectic::electric_flow(pic_domain, s.b, s.e, -1);
    s.b.data().array() += 0.13;
    symplectic::fill_periodic_ghosts(pic_domain, s.b);
    std::function<Eigen::VectorXd(const State&)> pack = [&](const State& input) {
        Eigen::VectorXd z(total);
        int r = 0;
        for (int i = 1; i <= 2; ++i)
            for (int j = 1; j <= 2; ++j)
                for (int k = 1; k <= 2; ++k)
                    for (int d = 0; d < 3; ++d, ++r) {
                        z[r] = input.e(i, j, k, d);
                        z[nf + r] = input.b(i, j, k, d);
                    }
        for (int d = 0; d < 3; ++d) {
            z[nx + d] = input.p.coord[d];
            z[nv + d] = input.p.velocity[d];
        }
        return z;
    };
    std::function<State(const Eigen::VectorXd&)> unpack = [&](const Eigen::VectorXd& z) {
        State result = state(domain);
        int r = 0;
        for (int i = 1; i <= 2; ++i)
            for (int j = 1; j <= 2; ++j)
                for (int k = 1; k <= 2; ++k)
                    for (int d = 0; d < 3; ++d, ++r) {
                        result.e(i, j, k, d) = z[r];
                        result.b(i, j, k, d) = z[nf + r];
                    }
        for (int d = 0; d < 3; ++d) {
            result.p.coord[d] = z[nx + d];
            result.p.velocity[d] = z[nv + d];
        }
        symplectic::fill_periodic_ghosts(pic_domain, result.e);
        symplectic::fill_periodic_ghosts(pic_domain, result.b);
        return result;
    };
    std::function<Eigen::MatrixXd(const State&)> bracket = [&](const State& input) {
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(total, total);
        for (int r = 0; r < nf; ++r) {
            Eigen::VectorXd basis = Eigen::VectorXd::Zero(total);
            basis[r] = 1;
            State unit = unpack(basis);
            Field3d magnetic = Field3d::Zero(domain.size(), 3);
            symplectic::electric_flow(pic_domain, magnetic, unit.e, 1);
            unit.b = magnetic;
            Eigen::VectorXd column = pack(unit);
            for (int t = 0; t < nf; ++t) {
                J(nf + t, r) = column[nf + t];
                J(r, nf + t) = -column[nf + t];
            }
            const Vector3R w = symplectic::gather<false>(pic_domain, unit.e, input.p.coord);
            for (int d = 0; d < 3; ++d) {
                J(nv + d, r) = -w[d];
                J(r, nv + d) = w[d];
            }
        }
        for (int d = 0; d < 3; ++d) {
            J(nx + d, nv + d) = 5;
            J(nv + d, nx + d) = -5;
        }
        const Vector3R b = symplectic::gather<true>(pic_domain, input.b, input.p.coord);
        for (int d = 0; d < 3; ++d) {
            const int a = (d + 1) % 3, c = (d + 2) % 3;
            J(nv + a, nv + c) = -5 * b[d];
            J(nv + c, nv + a) = 5 * b[d];
        }
        return J;
    };
    const Eigen::VectorXd z = pack(s);
    State advanced = s;
    step(domain, advanced, 0.04);
    Eigen::MatrixXd derivative(total, total);
    constexpr double epsilon = 1e-6;
    for (int r = 0; r < total; ++r) {
        Eigen::VectorXd plus = z, minus = z;
        plus[r] += epsilon;
        minus[r] -= epsilon;
        State sp = unpack(plus), sm = unpack(minus);
        step(domain, sp, 0.04);
        step(domain, sm, 0.04);
        derivative.col(r) = (pack(sp) - pack(sm)) / (2 * epsilon);
    }
    const Eigen::MatrixXd defect = derivative * bracket(s) * derivative.transpose() - bracket(advanced);
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

nlohmann::json settings() {
    return {{"Scheme", "symplectic"},
            {"Dt", 0.1},
            {"Dx", 1},
            {"Dy", 1},
            {"Dz", 1},
            {"NumCellsX", 3},
            {"NumCellsY", 3},
            {"NumCellsZ", 3},
            {"StartFromTime", 0},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

void parallel_deposits() {
    const int original_threads = omp_get_max_threads();
    std::function<std::unique_ptr<Simulation>(bool)> run = [&](bool parallel) {
        omp_set_num_threads(parallel ? 4 : 1);
        nlohmann::json cfg = settings();
        char name[] = "parallel_symplectic";
        char* argv[] = {name};
        std::unique_ptr<Simulation> base = build_simulation(
            cfg,
            {{"particles",
              nlohmann::json::array({check::species("Beam", -1, 1, 1, 512), check::species("Ions", 1, 25, 1, 512),
                                     check::species("Gas", 0, 2, 1, 512)})}},
            1, argv);
        base->init();
        SimulationSymplectic& sim = dynamic_cast<SimulationSymplectic&>(*base);
        std::mt19937_64 rng(71231);
        std::uniform_real_distribution<double> unit(0, 1);
        for (const Species::value_type& entry : sim.species) {
            const int count = entry.second->is_neutral() ? 5000 : 10000;
            for (int i = 0; i < count; ++i) {
                Vector3R x(3 * unit(rng), 3 * unit(rng), 3 * unit(rng));
                Vector3R v(0.2 * (unit(rng) - 0.5), 0.2 * (unit(rng) - 0.5), 0.2 * (unit(rng) - 0.5));
                if (i % 31 == 0) {
                    x[0] = 2.999;
                    v[0] = 0.1;
                }
                // Geometric stress case: a neutral traverses a complete
                // period in one step without generating electric current.
                if (entry.second->is_neutral() && i == 0) {
                    x = {0.5, 1.2, 0.4};
                    v = {31, 0, 0};
                }
                entry.second->add_particle(Particle(x, v));
            }
        }
        Field3d initial_density = Field3d::Zero(sim.domain.size(), 1);
        Field3d density = initial_density;
        Field3d current(sim.domain.size(), 3);
        sim.collect_charge_density(initial_density);
        Field3d old_density = initial_density;
        const Domain& pic_domain = sim.domain;
        for (int step = 1; step <= 3; ++step) {
            // Reuse scratch buffers after reducing the number of active
            // threads, then grow again. Old inactive deposits must not leak.
            omp_set_num_threads(parallel ? (step == 2 ? 2 : 4) : 1);
            sim.prepare_step(step);
            sim.make_step(step);
            sim.collect_charge_density(density);
            sim.collect_current(current);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k) {
                        const double rho = density(i + 1, j + 1, k + 1, 0);
                        check::near((rho - old_density(i + 1, j + 1, k + 1, 0)) / 0.1 +
                                        symplectic::divergence(pic_domain, current, i, j, k, false),
                                    0, 3e-13, "parallel local charge continuity");
                        check::near(symplectic::divergence(pic_domain, sim.electric, i, j, k, false) - rho,
                                    -initial_density(i + 1, j + 1, k + 1, 0), 3e-13, "parallel Gauss preservation");
                    }
            old_density = density;
        }
        sim.finalize();
        return base;
    };
    std::unique_ptr<Simulation> serial = run(false), parallel = run(true);
    const SimulationSymplectic& a = dynamic_cast<const SimulationSymplectic&>(*serial);
    const SimulationSymplectic& b = dynamic_cast<const SimulationSymplectic&>(*parallel);
    check::near((a.electric.data() - b.electric.data()).norm(), 0, 2e-13, "thread-count electric agreement");
    check::near((a.magnetic.data() - b.magnetic.data()).norm(), 0, 2e-13, "thread-count magnetic agreement");
    for (const Species::value_type& entry : a.species) {
        const ParticlesArray& sa = *entry.second;
        const ParticlesArray& sb = *b.species.at(entry.first);
        Field3d density_a(a.domain.size(), 1), density_b(b.domain.size(), 1);
        a.deposit_charge_density(sa, density_a);
        b.deposit_charge_density(sb, density_b);
        check::near((density_a.data() - density_b.data()).norm(), 0, 2e-13, "thread-count density agreement");
        check::near((sa.currentOnGrid.data() - sb.currentOnGrid.data()).norm(), 0, 2e-13,
                    "thread-count current agreement");
        for (int cell = 0; cell < sa.size(); ++cell) {
            check::require(sa(cell).size() == sb(cell).size(), "thread-count bucket population");
            for (size_t i = 0; i < sa(cell).size(); ++i) {
                check::vector(sa(cell)[i].coord, sb(cell)[i].coord, 3e-14, "thread-count position");
                check::vector(sa(cell)[i].velocity, sb(cell)[i].velocity, 3e-14, "thread-count velocity");
                if (sa.is_neutral() && sa(cell)[i].velocity[0] == 31)
                    check::vector(sa(cell)[i].coord, {0.8, 1.2, 0.4}, 2e-14, "neutral multi-period rebin");
            }
        }
    }
    std::unique_ptr<Simulation> repeated = run(true);
    const SimulationSymplectic& repeat = dynamic_cast<const SimulationSymplectic&>(*repeated);
    check::require((b.electric.data().array() == repeat.electric.data().array()).all(),
                   "fixed-thread electric result must be bitwise reproducible");
    check::require((b.magnetic.data().array() == repeat.magnetic.data().array()).all(),
                   "fixed-thread magnetic result must be bitwise reproducible");
    omp_set_num_threads(original_threads);
}

void simulation() {
    nlohmann::json cfg = settings();
    cfg["diagnostics"]["outputs"] = {{{"type", "energy_balance"}, {"interval", 1}},
                                     {{"type", "recovery"}, {"interval", 1}}};
    nlohmann::json particles = nlohmann::json{
        {"particles", nlohmann::json::array({check::species("Beam", -1, 1, 1, 1), check::species("Gas", 0, 2, 1, 1)})}};
    char name[] = "symplectic";
    char* argv[] = {name};
    std::unique_ptr<Simulation> base = build_simulation(cfg, particles, 1, argv);
    base->init();
    SimulationSymplectic& sim = dynamic_cast<SimulationSymplectic&>(*base);
    for (const Species::value_type& entry : sim.species)
        check::near(entry.second->currentOnGrid.norm(), 0, 0, "step-zero current must be initialized");
    ParticlesArray& beam = *sim.species.at("Beam");
    beam.add_particle(Particle({2.99, 0.8, 1.2}, {0.4, -0.1, 0.03}));
    ParticlesArray& gas = *sim.species.at("Gas");
    gas.add_particle(Particle({0.02, 1.2, 0.6}, {-1, 0.03, 0.02}));
    Field3d rho_old(sim.domain.size(), 1);
    sim.deposit_charge_density(beam, rho_old);
    sim.prepare_step(1);
    sim.make_step(1);
    const Domain& pic_domain = sim.domain;
    Field3d beam_density(sim.domain.size(), 1);
    sim.deposit_charge_density(beam, beam_density);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                const double drho = beam_density(i + 1, j + 1, k + 1, 0) - rho_old(i + 1, j + 1, k + 1, 0);
                check::near(drho / 0.1 + symplectic::divergence(pic_domain, beam.currentOnGrid, i, j, k, false), 0,
                            2e-14, "complete Simulation step charge continuity");
                check::near(symplectic::divergence(pic_domain, sim.electric, i, j, k, false) -
                                beam_density(i + 1, j + 1, k + 1, 0),
                            -rho_old(i + 1, j + 1, k + 1, 0), 3e-15, "Gauss defect preserved");
            }
    for (int cell = 0; cell < gas.size(); ++cell)
        for (const Particle& p : gas(cell)) {
            check::vector(p.coord, {2.92, 1.203, 0.602}, 1e-14, "neutral periodic drift");
            check::vector(p.velocity, {-1, 0.03, 0.02}, 0, "neutral velocity");
        }
    for (const Species::value_type& entry : sim.species)
        for (int cell = 0; cell < entry.second->size(); ++cell)
            for (const Particle& p : (*entry.second)(cell))
                check::require(sim.domain.sind(sim.domain.storage_cell_index(p.coord)[0],
                                               sim.domain.storage_cell_index(p.coord)[1],
                                               sim.domain.storage_cell_index(p.coord)[2]) == cell,
                               "particle not rebinned into its final cell");
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
    std::unique_ptr<Simulation> restart_base = build_simulation(restart_cfg, particles, 1, argv);
    restart_base->init();
    SimulationSymplectic& restart = dynamic_cast<SimulationSymplectic&>(*restart_base);
    restart.prepare_step(2);
    restart.make_step(2);
    check::near((restart.electric.data() - sim.electric.data()).norm(), 0, 1e-14, "restart electric field");
    check::near((restart.magnetic.data() - sim.magnetic.data()).norm(), 0, 1e-14, "restart magnetic field");
    for (const Species::value_type& entry : sim.species) {
        const ParticlesArray& restored = *restart.species.at(entry.first);
        for (int cell = 0; cell < entry.second->size(); ++cell) {
            check::require((*entry.second)(cell).size() == restored(cell).size(), "restart population");
            for (size_t p = 0; p < (*entry.second)(cell).size(); ++p) {
                check::vector(restored(cell)[p].coord, (*entry.second)(cell)[p].coord, 1e-14, "restart position");
                check::vector(restored(cell)[p].velocity, (*entry.second)(cell)[p].velocity, 1e-14, "restart velocity");
            }
        }
    }
    restart.make_diagnostic(0);
    restart.finalize();
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
            std::unique_ptr<Simulation> test = build_simulation(bad, particles, 1, argv);
            test->init();
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check::require(rejected, "unsupported symplectic configuration accepted");
    }
}
}   // namespace

int main() {
    return check::run([] {
        optimized_kernels();
        bucket_flows();
        periodic_ghosts();
        continuity();
        exact_subflows();
        poisson_map();
        convergence();
        plasma_stability_limit();
        simulation();
        parallel_deposits();
    });
}
