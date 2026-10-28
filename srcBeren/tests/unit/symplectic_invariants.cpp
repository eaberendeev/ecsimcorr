#include <filesystem>
#include <fstream>
#include <limits>
#include <numbers>

#include "simulation_symplectic.h"
#include "test_support.h"

using symplectic::WhitneyPIC;
namespace {
struct State {
    Field3d e, b;
    Particle p;
};
State state(const Domain& domain) {
    return {Field3d::Zero(domain.size(), 3), Field3d::Zero(domain.size(), 3),
            Particle(Vector3R(0.37, 0.41, 0.29), Vector3R(0.13, -0.07, 0.09))};
}
void step(const Domain& domain, State& s, double dt, double qm = -1, double qw = -0.2) {
    WhitneyPIC pic(domain);
    Field3d integral = Field3d::Zero(domain.size(), 3);
    auto drift = [&](int axis, double h) { pic.kinetic_flow(s.p, axis, h, qm, qw, s.e, s.b, integral); };
    drift(0, dt / 2);
    drift(1, dt / 2);
    drift(2, dt / 2);
    pic.magnetic_flow(s.e, s.b, dt / 2);
    pic.electric_flow(s.b, s.e, dt);
    s.p.velocity += dt * qm * pic.gather(s.e, s.p.coord, false);
    pic.magnetic_flow(s.e, s.b, dt / 2);
    drift(2, dt / 2);
    drift(1, dt / 2);
    drift(0, dt / 2);
    pic.sync(s.e);
    pic.sync(s.b);
}
void continuity() {
    Domain domain;
    domain.init({4, 3, 5}, {0.7, 1.1, 0.9});
    WhitneyPIC pic(domain);
    for (int axis = 0; axis < 3; ++axis)
        for (double velocity : {-19.3, -1.4, 0., 1.7, 21.2}) {
            auto s = state(domain);
            // Start exactly on a grid plane to exercise the negative path convention.
            s.p.coord = {0.7, 1.1, 0.9};
            s.p.velocity[axis] = velocity;
            const auto start = s.p.coord;
            Field3d old = Field3d::Zero(domain.size(), 1), rho = old;
            Field3d integral = Field3d::Zero(domain.size(), 3);
            pic.deposit_charge(old, start, -0.3);
            pic.kinetic_flow(s.p, axis, 1.3, -1, -0.3, s.e, s.b, integral);
            pic.deposit_charge(rho, s.p.coord, -0.3);
            const int g = domain.grid.ghost_cells();
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 5; ++k) {
                        const double change = rho(i + g, j + g, k + g, 0) - old(i + g, j + g, k + g, 0);
                        check::near(change + pic.divergence(integral, i, j, k, false), 0, 3e-14,
                                    "local continuity across cells and periodic wraps");
                    }
        }
}
void exact_subflows() {
    Domain domain;
    domain.init({4, 3, 5}, {0.7, 1.1, 0.9});
    WhitneyPIC pic(domain);
    auto s = state(domain);
    s.b.data().setConstant(0.6);
    const auto start = s.p;
    auto integral = Field3d::Zero(domain.size(), 3);
    pic.kinetic_flow(s.p, 0, 0.4, -2, 0, s.e, s.b, integral);
    check::vector(s.p.velocity,
                  start.velocity + Vector3R(0.4 * start.velocity[0], 0, 0).cross(Vector3R(0.6, 0.6, 0.6)) * -2.0, 1e-15,
                  "exact magnetic line integral and Lorentz sign");
    // Discrete curl and Whitney reconstruction commute inside a cell.
    for (int i = 0; i < s.e.capacity(); ++i) s.e[i] = std::sin(0.73 * i);
    pic.sync(s.e);
    auto curl = Field3d::Zero(domain.size(), 3);
    pic.electric_flow(curl, s.e, -1);
    const auto reconstructed = pic.gather(curl, start.coord, true);
    Vector3R finite_curl(0, 0, 0);
    for (int d = 0; d < 3; ++d) {
        const int a = (d + 1) % 3, b = (d + 2) % 3;
        auto plus = start.coord, minus = start.coord;
        plus[a] += 1e-6;
        minus[a] -= 1e-6;
        finite_curl[d] = (pic.gather(s.e, plus, false)[b] - pic.gather(s.e, minus, false)[b]) / 2e-6;
        plus = start.coord;
        minus = start.coord;
        plus[b] += 1e-6;
        minus[b] -= 1e-6;
        finite_curl[d] -= (pic.gather(s.e, plus, false)[a] - pic.gather(s.e, minus, false)[a]) / 2e-6;
    }
    check::vector(reconstructed, finite_curl, 5e-10, "Whitney curl compatibility");
    s.b = curl;
    pic.sync(s.b);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 5; ++k) check::near(pic.divergence(s.b, i, j, k, true), 0, 3e-15, "div curl E = 0");
    const auto original = s;
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
    WhitneyPIC pic(domain);
    auto s = state(domain);
    s.p.velocity = {0.1, 0, 0};
    double maximum = 0;
    for (int i = 0; i < steps; ++i) {
        step(domain, s, dt, -1, -1);
        const double energy = pic.field_energy(s.e) + 0.5 * s.p.velocity.dot(s.p.velocity);
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
double vacuum_error(double dt, int steps) {
    Domain domain;
    domain.init({8, 1, 1}, {0.7, 1.1, 0.9});
    auto s = state(domain);
    WhitneyPIC pic(domain);
    for (int i = 0; i < 8; ++i) s.e(i + 1, 1, 1, 1) = std::cos(2 * std::numbers::pi * i / 8);
    pic.sync(s.e);
    auto initial = s.e;
    auto curl = Field3d::Zero(domain.size(), 3);
    pic.electric_flow(curl, initial, -1);
    pic.sync(curl);
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
    WhitneyPIC pic(domain);
    constexpr int nf = 24, nx = 48, nv = 51, total = 54;
    auto s = state(domain);
    for (int i = 0; i < s.e.capacity(); ++i) s.e[i] = 0.03 * std::sin(0.71 * i);
    pic.sync(s.e);
    pic.electric_flow(s.b, s.e, -1);
    s.b.data().array() += 0.13;
    pic.sync(s.b);
    auto pack = [&](const State& input) {
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
    auto unpack = [&](const Eigen::VectorXd& z) {
        auto result = state(domain);
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
        pic.sync(result.e);
        pic.sync(result.b);
        return result;
    };
    auto bracket = [&](const State& input) {
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(total, total);
        for (int r = 0; r < nf; ++r) {
            Eigen::VectorXd basis = Eigen::VectorXd::Zero(total);
            basis[r] = 1;
            auto unit = unpack(basis);
            auto magnetic = Field3d::Zero(domain.size(), 3);
            pic.electric_flow(magnetic, unit.e, 1);
            unit.b = magnetic;
            auto column = pack(unit);
            for (int t = 0; t < nf; ++t) {
                J(nf + t, r) = column[nf + t];
                J(r, nf + t) = -column[nf + t];
            }
            const auto w = pic.gather(unit.e, input.p.coord, false);
            for (int d = 0; d < 3; ++d) {
                J(nv + d, r) = -w[d];
                J(r, nv + d) = w[d];
            }
        }
        for (int d = 0; d < 3; ++d) {
            J(nx + d, nv + d) = 5;
            J(nv + d, nx + d) = -5;
        }
        const auto b = pic.gather(input.b, input.p.coord, true);
        for (int d = 0; d < 3; ++d) {
            const int a = (d + 1) % 3, c = (d + 2) % 3;
            J(nv + a, nv + c) = -5 * b[d];
            J(nv + c, nv + a) = 5 * b[d];
        }
        return J;
    };
    const auto z = pack(s);
    auto advanced = s;
    step(domain, advanced, 0.04);
    Eigen::MatrixXd derivative(total, total);
    constexpr double epsilon = 1e-6;
    for (int r = 0; r < total; ++r) {
        auto plus = z, minus = z;
        plus[r] += epsilon;
        minus[r] -= epsilon;
        auto sp = unpack(plus), sm = unpack(minus);
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
            {"Collider", "None"},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}
void simulation() {
    auto cfg = settings();
    cfg["diagnostics"]["outputs"] = {{{"type", "energy_balance"}, {"interval", 1}},
                                     {{"type", "recovery"}, {"interval", 1}}};
    auto particles = nlohmann::json{
        {"particles", nlohmann::json::array({check::species("Beam", -1, 1, 1, 1), check::species("Gas", 0, 2, 1, 1)})}};
    char name[] = "symplectic";
    char* argv[] = {name};
    auto base = build_simulation(cfg, particles, 1, argv);
    base->init();
    auto& sim = dynamic_cast<SimulationSymplectic&>(*base);
    for (const auto& [name, sp] : sim.species)
        check::near(sp->currentOnGrid.norm(), 0, 0, "step-zero current must be initialized");
    auto& beam = *sim.species.at("Beam");
    beam.add_particle(Particle({2.99, 0.8, 1.2}, {0.4, -0.1, 0.03}));
    auto& gas = *sim.species.at("Gas");
    gas.add_particle(Particle({0.02, 1.2, 0.6}, {-1, 0.03, 0.02}));
    sim.deposit_density(beam);
    auto rho_old = beam.densityOnGrid;
    sim.prepare_step(1);
    sim.make_step(1);
    WhitneyPIC pic(sim.domain);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                const double drho = beam.densityOnGrid(i + 1, j + 1, k + 1, 0) - rho_old(i + 1, j + 1, k + 1, 0);
                check::near(drho / 0.1 + pic.divergence(beam.currentOnGrid, i, j, k, false), 0, 2e-14,
                            "complete Simulation step charge continuity");
                check::near(pic.divergence(sim.electric, i, j, k, false) - beam.densityOnGrid(i + 1, j + 1, k + 1, 0),
                            -rho_old(i + 1, j + 1, k + 1, 0), 3e-15, "Gauss defect preserved");
            }
    for (int cell = 0; cell < gas.size(); ++cell)
        for (const auto& p : gas.particlesData(cell)) {
            check::vector(p.coord, {2.92, 1.203, 0.602}, 1e-14, "neutral periodic drift");
            check::vector(p.velocity, {-1, 0.03, 0.02}, 0, "neutral velocity");
        }
    for (const auto& [key, sp] : sim.species)
        for (int cell = 0; cell < sp->size(); ++cell)
            for (const auto& p : sp->particlesData(cell))
                check::require(
                    sim.domain.sind(sim.domain.get_cell_index(p.coord)[0], sim.domain.get_cell_index(p.coord)[1],
                                    sim.domain.get_cell_index(p.coord)[2]) == cell,
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
    const auto original_path = std::filesystem::current_path();
    std::filesystem::create_directory("restart");
    std::filesystem::current_path("restart");
    auto restart_cfg = cfg;
    restart_cfg["StartFromTime"] = 1;
    restart_cfg["StartTimeStep"] = 1;
    restart_cfg["diagnostics"]["outputs"] = nlohmann::json::array();
    auto restart_base = build_simulation(restart_cfg, particles, 1, argv);
    restart_base->init();
    auto& restart = dynamic_cast<SimulationSymplectic&>(*restart_base);
    restart.prepare_step(2);
    restart.make_step(2);
    check::near((restart.electric.data() - sim.electric.data()).norm(), 0, 1e-14, "restart electric field");
    check::near((restart.magnetic.data() - sim.magnetic.data()).norm(), 0, 1e-14, "restart magnetic field");
    for (const auto& [name, sp] : sim.species) {
        const auto& restored = *restart.species.at(name);
        for (int cell = 0; cell < sp->size(); ++cell) {
            check::require(sp->particlesData(cell).size() == restored.particlesData(cell).size(), "restart population");
            for (size_t p = 0; p < sp->particlesData(cell).size(); ++p) {
                check::vector(restored.particlesData(cell)[p].coord, sp->particlesData(cell)[p].coord, 1e-14,
                              "restart position");
                check::vector(restored.particlesData(cell)[p].velocity, sp->particlesData(cell)[p].velocity, 1e-14,
                              "restart velocity");
            }
        }
    }
    restart.make_diagnostic(0);
    restart.finalize();
    std::filesystem::current_path(original_path);
    for (int variant = 0; variant < 8; ++variant) {
        auto bad = cfg;
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
            auto test = build_simulation(bad, particles, 1, argv);
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
        continuity();
        exact_subflows();
        poisson_map();
        convergence();
        simulation();
    });
}
