#include <Eigen/Dense>
#include <fstream>
#include <limits>
#include <random>

#include "DensityMeasurements.h"
#include "Diagnostic.h"
#include "ParticlesDiagnostic.h"
#include "simulation_ecsim_whitney.h"
#include "simulation_ecsim_whitney_frozen.h"
#include "sparse.h"
#include "test_support.h"
#include "whitney_shape.h"

namespace {
nlohmann::json settings(const std::string& scheme) {
    return {{"Scheme", scheme},
            {"Dt", 0.1},
            {"NumCellsX", 4},
            {"NumCellsY", 4},
            {"NumCellsZ", 4},
            {"Dx", 0.8},
            {"Dy", 1.1},
            {"Dz", 1.3},
            {"StartFromTime", 0},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

std::unique_ptr<Simulation> calculation(const std::string& scheme, bool plasma = true,
                                        const std::string& predictor = "", double dt = 0.1,
                                        const std::string& assembly = "reference") {
    char program[] = "ecsim_whitney";
    char* argv[] = {program};
    auto particles = nlohmann::json::array();
    if (plasma) {
        particles.push_back(check::species("electrons", -1, 1, 1, 8));
        particles.push_back(check::species("ions", 1, 25, 1, 8));
    }
    auto config = settings(scheme);
    config["Dt"] = dt;
    if (!predictor.empty())
        config["WhitneyOrbit"] = {{"predictor", predictor}};
    if (scheme == "ecsim_whitney_frozen")
        config["WhitneyOrbit"]["assembly"] = assembly;
    if (scheme == "ecsim_whitney_frozen_corr")
        config["WhitneyOrbit"] = {{"predictor", "ballistic"}, {"assembly", assembly}};
    if (scheme == "ecsim_whitney")
        config["WhitneyAssembly"] = assembly;
    auto sim = build_simulation(config, {{"particles", particles}}, 1, argv);
    sim->init();
    return sim;
}

void kernel_properties() {
    auto sim = calculation("ecsim_whitney", false);
    const auto& d = sim->domain;
    Field3d field(d.size(), 3);
    for (int i = 0; i < field.capacity(); ++i) field[i] = std::sin(0.31 * i);
    ecsimWhitney::sync_periodic(field, d);
    for (const Vector3R u : {Vector3R(1.2, 2.3, 0.7), Vector3R(3.9, 0.1, 3.8), Vector3R(0, 1, 2)}) {
        Vector3R position;
        for (int a = 0; a < 3; ++a) position[a] = u[a] * d.cell_size(a);
        Vector3R expected_e(0, 0, 0), expected_b(0, 0, 0);
        // Independent tensor-product definition, including exact grid planes.
        for (int c = 0; c < 3; ++c)
            for (int mask = 0; mask < 8; ++mask) {
                Vector3I node;
                double we = 1, wb = 1;
                for (int a = 0; a < 3; ++a) {
                    const int cell = static_cast<int>(std::floor(u[a])), bit = (mask >> a) & 1;
                    const double f = u[a] - cell;
                    node[a] = ecsimWhitney::periodic_index(cell + bit + d.grid.ghost_cells(), a, d);
                    we *= a == c ? (bit ? 0 : 1) : (bit ? f : 1 - f);
                    wb *= a != c ? (bit ? 0 : 1) : (bit ? f : 1 - f);
                }
                expected_e[c] += we * field(node[0], node[1], node[2], c);
                expected_b[c] += wb * field(node[0], node[1], node[2], c);
            }
        check::vector(ecsimWhitney::gather_electric(field, position, d), expected_e, 1e-13, "Whitney E tensor form");
        check::vector(ecsimWhitney::gather_magnetic(field, position, d), expected_b, 1e-13, "Whitney B tensor form");
        Field3d constant(d.size(), 3);
        constant.data().setConstant(1);
        check::vector(ecsimWhitney::gather_electric(constant, position, d), Vector3R(1), 1e-14,
                      "electric partition of unity");
        check::vector(ecsimWhitney::gather_magnetic(constant, position, d), Vector3R(1), 1e-14,
                      "magnetic partition of unity");
    }
    sim->finalize();
}

void compatible_gradient() {
    auto sim = calculation("ecsim_whitney", false);
    const auto& d = sim->domain;
    auto phi = [&](int i, int j, int k) {
        return std::sin(0.7 * (i % 4)) + std::cos(0.8 * (j % 4)) * std::sin(0.9 * (k % 4));
    };
    Field3d gradient(d.size(), 3);
    const int g = d.grid.ghost_cells();
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k) {
                gradient(i + g, j + g, k + g, 0) = (phi(i + 1, j, k) - phi(i, j, k)) / d.cell_size(0);
                gradient(i + g, j + g, k + g, 1) = (phi(i, j + 1, k) - phi(i, j, k)) / d.cell_size(1);
                gradient(i + g, j + g, k + g, 2) = (phi(i, j, k + 1) - phi(i, j, k)) / d.cell_size(2);
            }
    const Vector3R u(3.8, 1.2, 2.4);
    Vector3R position, exact(0, 0, 0);
    for (int a = 0; a < 3; ++a) position[a] = u[a] * d.cell_size(a);
    // Differentiate the nodal CIC scalar interpolant analytically. The edge
    // gather must commute with the gradient, including the periodic seam.
    for (int mask = 0; mask < 8; ++mask) {
        Vector3I node;
        for (int a = 0; a < 3; ++a) node[a] = static_cast<int>(std::floor(u[a])) + ((mask >> a) & 1);
        for (int c = 0; c < 3; ++c) {
            double derivative = ((mask >> c) & 1) ? 1 : -1;
            for (int a = 0; a < 3; ++a)
                if (a != c)
                    derivative *= ((mask >> a) & 1) ? u[a] - std::floor(u[a]) : 1 - u[a] + std::floor(u[a]);
            exact[c] += derivative * phi(node[0], node[1], node[2]) / d.cell_size(c);
        }
    }
    check::vector(ecsimWhitney::gather_electric(gradient, position, d), exact, 2e-14, "Whitney grad(CIC phi)");
    sim->finalize();
}

void orbit_continuity() {
    auto sim = calculation("ecsim_whitney", false);
    const auto& d = sim->domain;
    // Straight orbit crosses both an internal plane and the periodic seam.
    const Vector3R start(3.7, 1.8, 0.9), end(4.3, 2.4, 1.2);
    std::vector<double> cuts{0, 1};
    for (int a = 0; a < 3; ++a)
        for (int plane = static_cast<int>(std::floor(start[a])) + 1; plane < end[a]; ++plane)
            cuts.push_back((plane - start[a]) / (end[a] - start[a]));
    std::sort(cuts.begin(), cuts.end());
    Field3d current = Field3d::Zero(d.size(), 3), delta_charge = Field3d::Zero(d.size(), 1);
    auto charge = [&](const Vector3R& u, double sign) {
        ParticleShape<Shape, 2> shape;
        shape.fill_from_normalized(u);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k) {
                    const auto node = shape.start_ + Vector3I(i, j, k) + Vector3I(d.grid.ghost_cells());
                    delta_charge(ecsimWhitney::periodic_index(node[0], 0, d),
                                 ecsimWhitney::periodic_index(node[1], 1, d),
                                 ecsimWhitney::periodic_index(node[2], 2, d), 0) +=
                        sign * shape(i, 0) * shape(j, 1) * shape(k, 2);
                }
    };
    charge(start, -1);
    charge(end, 1);
    // Integral of W along each segment is quadratic: two Gauss points are
    // exact. dt=1; velocity is the physical endpoint displacement.
    for (size_t s = 1; s < cuts.size(); ++s)
        for (double gauss : {-1 / std::sqrt(3.0), 1 / std::sqrt(3.0)}) {
            const double half = 0.5 * (cuts[s] - cuts[s - 1]);
            const double t = 0.5 * (cuts[s] + cuts[s - 1]) + half * gauss;
            const auto u = start + t * (end - start);
            Vector3R position;
            for (int a = 0; a < 3; ++a) position[a] = u[a] * d.cell_size(a);
            for (const auto& w : ecsimWhitney::electric_weights(position, d))
                current[w.index] +=
                    half * w.weight * (end[w.component] - start[w.component]) * d.cell_size(w.component);
        }
    ecsimWhitney::sync_periodic(current, d);
    Vector3R physical_start, physical_end;
    for (int a = 0; a < 3; ++a) {
        physical_start[a] = start[a] * d.cell_size(a);
        physical_end[a] = end[a] * d.cell_size(a);
    }
    // Exercise production orbit weights against the independently assembled
    // current above, not just the quadrature written in this test.
    const auto orbit = ecsimWhitney::orbit_weights({physical_start, physical_start, physical_end, false}, d);
    Field3d production = Field3d::Zero(d.size(), 3);
    for (const auto& w : orbit.electric)
        production[w.index] += w.weight * (physical_end[w.component] - physical_start[w.component]);
    ecsimWhitney::sync_periodic(production, d);
    check::near((production - current).norm(), 0, 2e-14, "production orbit integral matches split quadrature");
    ecsimWhitney::sync_periodic(delta_charge, d);
    Field3d residual(d.size(), 1);
    spmv(sim->mesh.divE, current, residual);
    residual += delta_charge;
    check::near(residual.norm(), 0, 3e-14, "integrated Whitney current must satisfy CIC charge continuity");
    ecsimWhitney::OrbitWeights analytic;
    std::vector<double> scratch;
    ecsimWhitney::orbit_weights_cell({physical_start, physical_start, physical_end, false}, d, analytic, scratch);
    Field3d analytic_current = Field3d::Zero(d.size(), 3);
    for (const auto& w : analytic.electric)
        analytic_current[w.index] += w.weight * (physical_end[w.component] - physical_start[w.component]);
    ecsimWhitney::sync_periodic(analytic_current, d);
    spmv(sim->mesh.divE, analytic_current, residual);
    residual += delta_charge;
    check::near(residual.norm(), 0, 4e-14, "analytic orbit current satisfies independent CIC charge continuity");
    sim->finalize();
}

void two_leg_orbit_continuity() {
    auto sim = calculation("ecsim_whitney", false);
    const auto& d = sim->domain;
    // Two straight legs with different velocities, crossing an internal plane
    // and the periodic seam. The Boris turn is at the middle point.
    const Vector3R start(3.7, 1.8, 0.9), middle(4.0, 2.1, 1.05), end(4.3, 2.4, 1.2);
    const double dt = 1.0;
    auto physical = [&](const Vector3R& u) {
        Vector3R result;
        for (int a = 0; a < 3; ++a) result[a] = u[a] * d.cell_size(a);
        return result;
    };
    Field3d current = Field3d::Zero(d.size(), 3), delta_charge = Field3d::Zero(d.size(), 1);
    auto charge = [&](const Vector3R& u, double sign) {
        ParticleShape<Shape, 2> shape;
        shape.fill_from_normalized(u);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k) {
                    const auto node = shape.start_ + Vector3I(i, j, k) + Vector3I(d.grid.ghost_cells());
                    delta_charge(ecsimWhitney::periodic_index(node[0], 0, d),
                                 ecsimWhitney::periodic_index(node[1], 1, d),
                                 ecsimWhitney::periodic_index(node[2], 2, d), 0) +=
                        sign * shape(i, 0) * shape(j, 1) * shape(k, 2);
                }
    };
    charge(start, -1);
    charge(end, 1);
    const Vector3R pstart = physical(start), pmiddle = physical(middle), pend = physical(end);
    const Vector3R v1 = (pmiddle - pstart) / (0.5 * dt), v2 = (pend - pmiddle) / (0.5 * dt);
    ecsimWhitney::deposit_orbit_current(current, 1.0, {{pstart, pmiddle, v1, 0.5}, {pmiddle, pend, v2, 0.5}}, d);
    ecsimWhitney::sync_periodic(current, d);
    ecsimWhitney::sync_periodic(delta_charge, d);
    Field3d residual(d.size(), 1);
    spmv(sim->mesh.divE, current, residual);
    residual += delta_charge * (1 / dt);
    check::near(residual.norm(), 0, 3e-14, "two-leg Whitney current must satisfy CIC charge continuity");
    sim->finalize();
}

void analytic_orbit_properties() {
    using namespace ecsimWhitney;
    auto sim = calculation("ecsim_whitney_frozen", false);
    const auto& d = sim->domain;
    OrbitWeights optimized;
    std::vector<double> cuts;
    std::mt19937 random(8231);
    std::uniform_real_distribution<double> start_distribution(-4, 8), shift_distribution(-3, 3);
    auto world = [&](const Vector3R& u) {
        Vector3R position;
        for (int a = 0; a < 3; ++a) position[a] = u[a] * d.cell_size(a);
        return position;
    };
    auto compare = [&](const FrozenPath& path) {
        const auto reference = orbit_weights(path, d);
        orbit_weights_cell(path, d, optimized, cuts);
        for (bool electric : {true, false}) {
            const auto& actual = electric ? optimized.electric : optimized.magnetic;
            const auto& expected = electric ? reference.electric : reference.magnetic;
            std::map<int, double> difference;
            Vector3R sums(0, 0, 0);
            for (const auto& w : actual) {
                check::require(std::isfinite(w.weight) && w.weight >= -1e-14, "finite nonnegative orbit weights");
                difference[w.index] += w.weight;
                sums[w.component] += w.weight;
            }
            for (const auto& w : expected) difference[w.index] -= w.weight;
            for (const auto& [index, value] : difference)
                check::near(value, 0, 5e-14, "analytic weights equal independent split Gauss quadrature");
            check::vector(sums, Vector3R(1), 5e-14, "analytic orbit partition of unity");
        }
    };
    // Zero length, exact planes, simultaneous corner crossings, reverse
    // motion, negative coordinates, and several unwrapped periodic turns.
    for (const Vector3R end :
         {Vector3R(0, 1, 2), Vector3R(1, 1, 2), Vector3R(4, 5, 6), Vector3R(-4, -3, -2), Vector3R(12.2, -8.4, 5.1)}) {
        const auto start = world(Vector3R(0, 1, 2)), last = world(end);
        compare({start, start, last, false});
        compare({last, last, start, false});
    }
    for (int sample = 0; sample < 100; ++sample) {
        Vector3R a, b, c;
        for (int axis = 0; axis < 3; ++axis) {
            a[axis] = start_distribution(random);
            b[axis] = a[axis] + shift_distribution(random);
            c[axis] = b[axis] + shift_distribution(random);
        }
        compare({world(a), world(b), world(c), sample % 2 == 0});
    }
    for (const double bad : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN(), 1e10}) {
        bool rejected = false;
        try {
            orbit_weights_cell({Vector3R(0.0), Vector3R(0.0), Vector3R(bad, 0, 0), false}, d, optimized, cuts);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check::require(rejected, "analytic orbit rejects invalid coordinates before integer conversion");
    }
    sim->finalize();
}

void frozen_orbit_limits_and_predictors() {
    using namespace ecsimWhitney;
    auto sim = calculation("ecsim_whitney_frozen", false);
    const auto& d = sim->domain;
    Field3d field(d.size(), 3), zero = Field3d::Zero(d.size(), 3), constant(d.size(), 3);
    for (int i = 0; i < field.capacity(); ++i) field[i] = std::sin(0.17 * i);
    sync_periodic(field, d);
    constant.data().setConstant(0.3);
    const Vector3R start(0.9, 1.4, 1.6), end(1.3, 1.7, 2.1), velocity(0.2, -0.1, 0.05);
    const auto midpoint = 0.5 * (start + end);
    const auto weights = orbit_weights({start, midpoint, end, false}, d);
    // All points are in one cell. Independent Simpson integration is exact
    // for these quadratic E and linear B interpolants.
    check::vector(gather_weights(field, weights.electric),
                  (gather_electric(field, start, d) + 4.0 * gather_electric(field, midpoint, d) +
                   gather_electric(field, end, d)) /
                      6.0,
                  2e-14, "orbit E equals exact Simpson average inside a cell");
    check::vector(gather_weights(field, weights.magnetic),
                  (gather_magnetic(field, start, d) + gather_magnetic(field, end, d)) / 2.0, 2e-14,
                  "orbit B equals exact linear average inside a cell");
    for (const Vector3R last : {start, Vector3R(-1.1, 3.8, 5.5)}) {
        const auto forward = orbit_weights({start, start, last, false}, d);
        const auto reverse = orbit_weights({last, last, start, false}, d);
        check::vector(gather_weights(field, forward.electric), gather_weights(field, reverse.electric), 3e-14,
                      "reversed multi-cell orbit has the same averaged E gather");
        check::vector(gather_weights(field, forward.magnetic), gather_weights(field, reverse.magnetic), 3e-14,
                      "reversed multi-cell orbit has the same averaged B gather");
        check::vector(gather_weights(constant, forward.electric), Vector3R(0.3), 2e-14,
                      "orbit electric partition of unity");
        check::vector(gather_weights(constant, forward.magnetic), Vector3R(0.3), 2e-14,
                      "orbit magnetic partition of unity");
        if ((last - start).norm() == 0) {
            check::vector(gather_weights(field, forward.electric), gather_electric(field, start, d), 1e-14,
                          "zero-length orbit reduces to the point electric gather");
            check::vector(gather_weights(field, forward.magnetic), gather_magnetic(field, start, d), 1e-14,
                          "zero-length orbit reduces to the point magnetic gather");
        }
    }
    const double dt = 0.2, qm = -0.7;
    for (const auto predictor : {OrbitPredictor::Ballistic, OrbitPredictor::OldFields, OrbitPredictor::HalfDrift,
                                 OrbitPredictor::CenteredFields}) {
        const auto path = predict_path(start, velocity, qm, dt, constant, zero, d, predictor);
        const bool centered = predictor == OrbitPredictor::CenteredFields;
        const auto acceleration =
            predictor == OrbitPredictor::Ballistic || centered ? Vector3R(0.0) : Vector3R(0.3 * qm);
        check::vector(path.end, start + dt * velocity + (0.5 * dt * dt) * acceleration, 1e-14,
                      "constant-E predictor has the analytic accelerated endpoint");
        check::require(path.split == (predictor == OrbitPredictor::HalfDrift || centered),
                       "wrong frozen path topology");
        if (centered)
            check::vector(path.middle, start + (0.5 * dt) * velocity - (dt * dt / 8) * Vector3R(0.3 * qm), 1e-14,
                          "centered predictor has the analytic constant-E half-time position");
        const auto free_path = predict_path(start, velocity, qm, dt, zero, zero, d, predictor);
        const auto free_weights = orbit_weights(free_path, d);
        const auto ballistic = orbit_weights({start, start, start + dt * velocity, false}, d);
        check::vector(gather_weights(field, free_weights.electric), gather_weights(field, ballistic.electric), 2e-14,
                      "all predictors coincide for free streaming");
    }
    sim->finalize();
}

void seed_plasma(Simulation& sim) {
    std::mt19937 random(731);
    std::uniform_real_distribution<double> uniform(0, 1);
    std::normal_distribution<double> gaussian(0, 0.12);
    for (auto& [name, sp] : sim.species)
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    for (int p = 0; p < 8; ++p) {
                        const double x = uniform(random), y = uniform(random), z = uniform(random);
                        const double vx = gaussian(random), vy = gaussian(random), vz = gaussian(random);
                        sp->add_particle(
                            Particle(Vector3R((i + x) * sim.domain.cell_size(0), (j + y) * sim.domain.cell_size(1),
                                              (k + z) * sim.domain.cell_size(2)),
                                     Vector3R(vx, vy, vz) / std::sqrt(sp->mass())));
                    }
}

void response_matches_boris_work(const std::string& predictor = "", const std::string& assembly = "reference") {
    auto sim =
        calculation(predictor.empty() ? "ecsim_whitney" : "ecsim_whitney_frozen", true, predictor, 0.1, assembly);
    seed_plasma(*sim);
    auto& s = dynamic_cast<SimulationEcsimWhitney&>(*sim);
    // Nonzero B exercises off-diagonal alpha signs, not just its diagonal.
    for (int i = 0; i < s.fieldBn.capacity(); ++i) s.fieldBn[i] = 0.1 * std::sin(0.07 * i);
    ecsimWhitney::sync_periodic(s.fieldBn, s.domain);
    s.prepare_step(1);
    s.drift_particles_and_assemble_response();
    for (int i = 0; i < s.fieldEp.capacity(); ++i) {
        s.fieldEp[i] = 0.2 * std::cos(0.11 * i);
        s.fieldE_external[i] = 0.05 * std::sin(0.13 * i);
    }
    ecsimWhitney::sync_periodic(s.fieldEp, s.domain);
    ecsimWhitney::sync_periodic(s.fieldE_external, s.domain);
    double before = 0;
    for (const auto& [name, sp] : s.species) before += sp->get_kinetic_energy();
    s.advance_particle_velocities();
    const Field3d electric = s.fieldEp + s.fieldE_external;
    const Field3d predicted = s.fieldJp + (2 / s.dt_) * (s.mesh.Lmat2 * electric);
    Field3d deposited(s.domain.size(), 3);
    s.collect_current(deposited);
    check::near((predicted - deposited).norm(), 0, 2e-12, "L/Jp must reproduce Boris midpoint current and aliases");
    double after = 0;
    for (const auto& [name, sp] : s.species) after += sp->get_kinetic_energy();
    const double work = s.dt_ * dot_product_sum(electric, deposited, s.bc_handler.active_range(s.domain.grid));
    check::near(after - before, work, 2e-12, "particle kinetic change equals adjoint deposited electric work");
    s.finalize();
}

void cell_assembly_matches_reference(const std::string& predictor) {
    const std::string scheme = predictor.empty() ? "ecsim_whitney" : "ecsim_whitney_frozen";
    auto reference = calculation(scheme, true, predictor, 0.1, "reference");
    auto optimized = calculation(scheme, true, predictor, 0.1, "cell");
    for (auto* sim : {reference.get(), optimized.get()}) {
        seed_plasma(*sim);
        auto& s = dynamic_cast<SimulationEcsim&>(*sim);
        for (int i = 0; i < s.fieldBn.capacity(); ++i) {
            s.fieldBn[i] = 0.1 * std::sin(0.07 * i);
            s.fieldEn[i] = 0.1 * std::cos(0.11 * i);
        }
        ecsimWhitney::sync_periodic(s.fieldBn, s.domain);
        ecsimWhitney::sync_periodic(s.fieldEn, s.domain);
        s.prepare_step(1);
        s.drift_particles_and_assemble_response();
    }
    auto& a = dynamic_cast<SimulationEcsim&>(*reference);
    auto& b = dynamic_cast<SimulationEcsim&>(*optimized);
    check::near((a.mesh.Lmat2 - b.mesh.Lmat2).norm(), 0, 2e-14, "cell sums preserve the complete orbit L matrix");
    check::near((a.fieldJp - b.fieldJp).norm(), 0, 2e-14, "cell assembly preserves the orbit predictor current");
    // Compare the actual push as well: matching L alone would miss races in
    // the velocity update or its adjoint current scatter.
    for (auto* s : {&a, &b}) {
        for (int i = 0; i < s->fieldEp.capacity(); ++i) {
            s->fieldEp[i] = 0.2 * std::cos(0.11 * i);
            s->fieldE_external[i] = 0.05 * std::sin(0.13 * i);
        }
        ecsimWhitney::sync_periodic(s->fieldEp, s->domain);
        ecsimWhitney::sync_periodic(s->fieldE_external, s->domain);
        s->advance_particle_velocities();
    }
    for (const auto& [name, sp] : a.species) {
        const auto& other = *b.species.at(name);
        check::near((sp->currentOnGrid - other.currentOnGrid).norm(), 0, 2e-13,
                    "parallel push preserves the complete adjoint current");
        for (int cell = 0; cell < sp->size(); ++cell) {
            const auto& particles = sp->particlesData(cell);
            const auto& expected = other.particlesData(cell);
            check::require(particles.size() == expected.size(), "parallel push preserves particle buckets");
            for (size_t p = 0; p < particles.size(); ++p)
                check::vector(particles[p].velocity, expected[p].velocity, 2e-14,
                              "parallel and reference Boris velocities agree");
        }
    }
    a.finalize();
    b.finalize();
}

void cell_predictor_failure_propagates() {
    auto sim = calculation("ecsim_whitney_frozen", true, "old_fields", 0.1, "cell");
    seed_plasma(*sim);
    auto& s = dynamic_cast<SimulationEcsim&>(*sim);
    s.fieldE_external.data().setConstant(std::numeric_limits<double>::infinity());
    s.prepare_step(1);
    bool rejected = false;
    try {
        s.drift_particles_and_assemble_response();
    } catch (const std::runtime_error& error) {
        rejected = std::string(error.what()).find("non-finite") != std::string::npos;
    }
    check::require(rejected, "parallel orbit failure must propagate to the simulation caller");
    s.finalize();
}

void cached_response_lifetime() {
    auto sim = calculation("ecsim_whitney_frozen", true, "ballistic", 0.1, "cell");
    seed_plasma(*sim);
    auto& s = dynamic_cast<SimulationEcsim&>(*sim);
    s.prepare_step(1);
    s.drift_particles_and_assemble_response();
    s.advance_particle_velocities();
    bool rejected = false;
    try {
        s.advance_particle_velocities();
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    check::require(rejected, "cached orbit allocations must not make a consumed response valid again");
    // Injection changes the response cardinality. The next assembly must
    // overwrite cached weights in the new bucket/particle traversal order.
    for (auto& [name, sp] : s.species) sp->add_particle(Particle(Vector3R(0.7, 1.2, 1.6), Vector3R(0.02, -0.03, 0.01)));
    s.prepare_step(2);
    s.drift_particles_and_assemble_response();
    for (int i = 0; i < s.fieldEp.capacity(); ++i) s.fieldEp[i] = 0.1 * std::cos(0.11 * i);
    ecsimWhitney::sync_periodic(s.fieldEp, s.domain);
    double before = 0;
    for (const auto& [name, sp] : s.species) before += sp->get_kinetic_energy();
    s.advance_particle_velocities();
    Field3d deposited(s.domain.size(), 3);
    s.collect_current(deposited);
    const Field3d predicted = s.fieldJp + (2 / s.dt_) * (s.mesh.Lmat2 * s.fieldEp);
    check::near((predicted - deposited).norm(), 0, 2e-12, "rebuilt cached response reproduces the new Boris current");
    double after = 0;
    for (const auto& [name, sp] : s.species) after += sp->get_kinetic_energy();
    const double work = s.dt_ * dot_product_sum(s.fieldEp, deposited, s.bc_handler.active_range(s.domain.grid));
    check::near(after - before, work, 2e-12, "cached response after injection preserves electric work");
    s.finalize();
}

void initialize_gauss_field(SimulationEcsim& s) {
    const auto rho = measure_charge_density(s.species, s.domain, s.bc_handler, ShapeType::Linear);
    const int g = s.domain.grid.ghost_cells();
    auto index = [](int i, int j, int k) { return 16 * ((i + 4) % 4) + 4 * ((j + 4) % 4) + (k + 4) % 4; };
    Eigen::MatrixXd laplacian = Eigen::MatrixXd::Zero(64, 64);
    Eigen::VectorXd rhs(64);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k) {
                const int row = index(i, j, k);
                rhs[row] = rho(i + g, j + g, k + g, 0);
                for (int a = 0; a < 3; ++a) {
                    const double h = s.domain.cell_size(a);
                    Vector3I lo(i, j, k), hi(i, j, k);
                    --lo[a];
                    ++hi[a];
                    laplacian(row, row) -= 2 / (h * h);
                    laplacian(row, index(lo[0], lo[1], lo[2])) += 1 / (h * h);
                    laplacian(row, index(hi[0], hi[1], hi[2])) += 1 / (h * h);
                }
            }
    check::near(rhs.sum(), 0, 2e-13, "periodic Poisson requires zero net charge");
    // Fix the potential gauge. E=grad(phi), div(E)=rho with this sign.
    laplacian.row(0).setZero();
    laplacian(0, 0) = 1;
    rhs[0] = 0;
    const Eigen::VectorXd potential = laplacian.colPivHouseholderQr().solve(rhs);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                for (int a = 0; a < 3; ++a) {
                    Vector3I hi(i, j, k);
                    ++hi[a];
                    s.fieldEn(i + g, j + g, k + g, a) =
                        (potential[index(hi[0], hi[1], hi[2])] - potential[index(i, j, k)]) / s.domain.cell_size(a);
                }
    ecsimWhitney::sync_periodic(s.fieldEn, s.domain);
    Field3d divergence(s.domain.size(), 1);
    spmv(s.mesh.divE, s.fieldEn, divergence);
    check::near((divergence - rho).norm(), 0, 5e-13, "initial electric field must satisfy discrete Gauss law");
}

void corr_charge_continuity(const std::string& scheme) {
    auto sim = calculation(scheme, true);
    seed_plasma(*sim);
    auto& s = dynamic_cast<SimulationEcsim&>(*sim);
    initialize_gauss_field(s);
    s.prepare_step(1);
    const auto before = measure_charge_density(s.species, s.domain, s.bc_handler, ShapeType::Linear);
    s.make_step(1);
    const auto after = measure_charge_density(s.species, s.domain, s.bc_handler, ShapeType::Linear);
    Field3d residual(s.domain.size(), 1), divergence(s.domain.size(), 1);
    spmv(s.mesh.divE, s.fieldJe, divergence);
    residual = divergence + (after - before) * (1 / s.dt_);
    check::near(residual.norm(), 0, 1e-10, scheme + " path current must satisfy CIC charge continuity");
    Field3d gauss(s.domain.size(), 1);
    spmv(s.mesh.divE, s.fieldEn, gauss);
    gauss -= after;
    check::near(gauss.norm(), 0, 1e-9, scheme + " corrected field must preserve discrete Gauss law");
    check::require(s.fieldEn.data().allFinite() && s.fieldBn.data().allFinite(),
                   scheme + " produced non-finite fields");
    s.finalize();
}

void corr_energy_conservation(const std::string& scheme) {
    auto sim = calculation(scheme, true);
    seed_plasma(*sim);
    auto& s = dynamic_cast<SimulationEcsim&>(*sim);
    initialize_gauss_field(s);
    const auto range = s.bc_handler.active_range(s.domain.grid);
    auto total_energy = [&] {
        double result =
            0.5 * dot_product_sum(s.fieldEn, s.fieldEn, range) + 0.5 * dot_product_sum(s.fieldBn, s.fieldBn, range);
        for (const auto& [name, sp] : s.species) result += sp->get_kinetic_energy();
        return result;
    };
    double initial_energy = 0;
    for (int step = 1; step <= 20; ++step) {
        s.prepare_step(step);
        if (step == 1)
            initial_energy = total_energy();
        s.make_step(step);
        check::require(s.fieldEn.data().allFinite() && s.fieldBn.data().allFinite(),
                       scheme + " produced non-finite fields during energy test");
    }
    const double final_energy = total_energy();
    check::near(final_energy, initial_energy, 1e-6 * std::max(1.0, std::abs(initial_energy)),
                scheme + " must conserve total energy over 20 steps");
    s.finalize();
}

void compare_schemes() {
    std::ofstream csv("ecsim_whitney_comparison.csv");
    csv << "case,scheme,predictor,step,electric_rms,gauss_rms,gauss_drift_rms,continuity_rms,energy_error\n";
    const std::vector<std::pair<std::string, std::string>> variants = {{"ecsim", ""},
                                                                       {"ecsim_whitney", ""},
                                                                       {"ecsim_whitney_frozen", "ballistic"},
                                                                       {"ecsim_whitney_frozen", "old_fields"},
                                                                       {"ecsim_whitney_frozen", "half_drift"},
                                                                       {"ecsim_whitney_frozen", "centered_fields"}};
    struct ComparisonCase {
        std::string name;
        double dt, bz;
    };
    const std::vector<ComparisonCase> cases = {{"dt01", 0.1, 0}, {"dt04", 0.4, 0}, {"magnetized_dt04", 0.4, 2}};
    for (const auto& experiment : cases)
        for (const auto& [scheme, predictor] : variants) {
            auto sim = calculation(scheme, true, predictor, experiment.dt);
            seed_plasma(*sim);
            auto& s = dynamic_cast<SimulationEcsim&>(*sim);
            initialize_gauss_field(s);
            for (int i = 0; i < s.fieldBn.capacity(); ++i) s.fieldBn[i] = i % 3 == 2 ? experiment.bz : 0;
            const auto range = s.bc_handler.active_range(s.domain.grid);
            Field3d initial_defect(s.domain.size(), 1);
            spmv(s.mesh.divE, s.fieldEn, initial_defect);
            initial_defect -= measure_charge_density(s.species, s.domain, s.bc_handler, ShapeType::Linear);
            double initial_energy = 0;
            for (int i = range.start[0]; i < range.end[0]; ++i)
                for (int j = range.start[1]; j < range.end[1]; ++j)
                    for (int k = range.start[2]; k < range.end[2]; ++k)
                        for (int c = 0; c < 3; ++c)
                            initial_energy +=
                                0.5 * (std::pow(s.fieldEn(i, j, k, c), 2) + std::pow(s.fieldBn(i, j, k, c), 2));
            for (const auto& [name, sp] : s.species) initial_energy += sp->get_kinetic_energy();
            for (int step = 1; step <= 20; ++step) {
                const auto old_charge = measure_charge_density(s.species, s.domain, s.bc_handler, ShapeType::Linear);
                s.prepare_step(step);
                s.make_step(step);
                check::require(s.fieldEn.data().allFinite() && s.fieldBn.data().allFinite(),
                               "non-finite plasma fields");
                Field3d divergence(s.domain.size(), 1);
                spmv(s.mesh.divE, s.fieldEn, divergence);
                const auto charge = measure_charge_density(s.species, s.domain, s.bc_handler, ShapeType::Linear);
                // Use the current that actually drives Maxwell. This also defines
                // a fair baseline for ordinary ECSIM, whose currentOnGrid is only
                // populated later by its separate diagnostic path.
                const Field3d response_current = s.fieldJp + (2 / s.dt_) * (s.mesh.Lmat2 * s.fieldEp);
                Field3d continuity(s.domain.size(), 1);
                spmv(s.mesh.divE, response_current, continuity);
                continuity += (charge - old_charge) * (1 / s.dt_);
                double gauss = 0, drift = 0, e2 = 0, energy = 0, continuity2 = 0;
                for (int i = range.start[0]; i < range.end[0]; ++i)
                    for (int j = range.start[1]; j < range.end[1]; ++j)
                        for (int k = range.start[2]; k < range.end[2]; ++k) {
                            const double defect = divergence(i, j, k, 0) - charge(i, j, k, 0);
                            gauss += defect * defect;
                            drift += std::pow(defect - initial_defect(i, j, k, 0), 2);
                            continuity2 += std::pow(continuity(i, j, k, 0), 2);
                            for (int c = 0; c < 3; ++c) {
                                e2 += std::pow(s.fieldEn(i, j, k, c), 2);
                                energy +=
                                    0.5 * (std::pow(s.fieldEn(i, j, k, c), 2) + std::pow(s.fieldBn(i, j, k, c), 2));
                            }
                        }
                for (const auto& [name, sp] : s.species) energy += sp->get_kinetic_energy();
                check::require(std::isfinite(energy) && std::isfinite(gauss), "non-finite plasma metrics");
                if (scheme != "ecsim")
                    check::near(energy, initial_energy, 2e-8, "periodic Whitney ECSIM must conserve total energy");
                csv << experiment.name << ',' << scheme << ',' << predictor << ',' << step << ',' << std::sqrt(e2 / 192)
                    << ',' << std::sqrt(gauss / 64) << ',' << std::sqrt(drift / 64) << ','
                    << std::sqrt(continuity2 / 64) << ',' << energy - initial_energy << '\n';
            }
            s.finalize();
        }
    check::require(csv.good(), "comparison CSV write failed");
}

void ballistic_wall_path_clipping() {
    using namespace ecsimWhitney;
    Domain d;
    d.init(Vector3I(4, 4, 4), Vector3R(1.0, 1.0, 1.0));
    d.periodic = {false, false, false};

    PathPolicy absorbing;
    absorbing.enable(Face::XMIN);
    absorbing.classify = [](const Particle&, Face face) { return ParticleFateResult{ParticleFate::Removed, face}; };

    const Vector3R start(0.5, 2.0, 2.0), velocity(-10.0, 0.0, 0.0);
    const FrozenPath free{start, start, start + 0.2 * velocity, false};
    Particle particle(free.end, velocity);
    ParticleFateResult fate;
    const FrozenPath clipped = clip_ballistic_path(free, particle, absorbing, d, fate);
    check::require(fate.fate == ParticleFate::Removed, "wall path records absorption");
    check::require(clipped.wall == Face::XMIN, "wall path records the first face");
    check::near(clipped.duration, 0.25, 0.0, "wall path covers exactly the incident time fraction");
    check::near(clipped.end.x(), 0.0, 0.0, "wall path ends on the wall plane");
    check::near(clipped.end.y(), start.y(), 0.0, "wall path keeps the transverse coordinate");
    check::require(!clipped.split, "absorbed path has one chord");

    const auto weights = orbit_weights(clipped, d);
    Vector3R electric_sum(0, 0, 0), magnetic_sum(0, 0, 0);
    for (const auto& w : weights.electric) electric_sum[w.component] += w.weight;
    for (const auto& w : weights.magnetic) magnetic_sum[w.component] += w.weight;
    check::vector(electric_sum, Vector3R(0.25), 2e-14, "absorbed orbit E weights sum to t_hit");
    check::vector(magnetic_sum, Vector3R(0.25), 2e-14, "absorbed orbit B weights sum to t_hit");

    // A particle entering the box from the wall must not be clipped.
    const FrozenPath inward{Vector3R(0.0, 2.0, 2.0), Vector3R(0.0, 2.0, 2.0), Vector3R(0.2, 2.0, 2.0), false};
    const FrozenPath untouched =
        clip_ballistic_path(inward, Particle(inward.end, Vector3R(2, 0, 0)), absorbing, d, fate);
    check::require(fate.fate == ParticleFate::Unhandled, "inward path has no wall event");
    check::near(untouched.duration, 1.0, 0.0, "inward path covers the full step");

    // Reflected paths are supported by the kernel but rejected by the
    // simulation because no wall surface charge is modeled here.
    PathPolicy reflecting = absorbing;
    reflecting.classify = [](const Particle&, Face face) { return ParticleFateResult{ParticleFate::Reflected, face}; };
    fate = {};
    const FrozenPath reflected = clip_ballistic_path(free, particle, reflecting, d, fate);
    check::require(fate.fate == ParticleFate::Reflected, "reflecting wall records reflection");
    check::require(reflected.split, "reflected path has two legs");
    check::near(reflected.first_fraction, 0.25, 0.0, "reflected path uses t_hit for the incident leg");
    check::near(reflected.end.x(), 1.5, 1e-15, "reflected endpoint mirrors the free endpoint");
    const auto reflected_weights = orbit_weights(reflected, d);
    electric_sum = Vector3R(0, 0, 0);
    magnetic_sum = Vector3R(0, 0, 0);
    for (const auto& w : reflected_weights.electric) electric_sum[w.component] += w.weight;
    for (const auto& w : reflected_weights.magnetic) magnetic_sum[w.component] += w.weight;
    check::vector(electric_sum, Vector3R(1), 3e-14, "reflected orbit E weights cover the full step");
    check::vector(magnetic_sum, Vector3R(1), 3e-14, "reflected orbit B weights cover the full step");

    // A path ending exactly on the wall is clipped at t=1 and absorbed;
    // the path still contributes its full-step current.
    absorbing.enable(Face::XMAX);
    const FrozenPath touching{Vector3R(3.5, 2, 2), Vector3R(3.5, 2, 2), Vector3R(4, 2, 2), false};
    fate = {};
    const FrozenPath touched =
        clip_ballistic_path(touching, Particle(touching.end, Vector3R(5, 0, 0)), absorbing, d, fate);
    check::require(fate.fate == ParticleFate::Removed, "path reaching the wall at the endpoint is absorbed");
    check::near(touched.duration, 1.0, 0.0, "endpoint wall contact keeps the full-step duration");

    // The cell-local accumulator must reproduce the reference split quadrature
    // for the clipped and reflected topologies as well.
    OrbitWeights optimized;
    std::vector<double> cuts;
    const auto compare_weights = [&](const FrozenPath& path, const std::string& label) {
        const auto reference = orbit_weights(path, d);
        orbit_weights_cell(path, d, optimized, cuts);
        for (bool electric : {true, false}) {
            const auto& actual = electric ? optimized.electric : optimized.magnetic;
            const auto& expected = electric ? reference.electric : reference.magnetic;
            std::map<int, double> difference;
            for (const auto& w : actual) difference[w.index] += w.weight;
            for (const auto& w : expected) difference[w.index] -= w.weight;
            for (const auto& [index, value] : difference)
                check::near(value, 0, 5e-14, label + " cell weights equal reference quadrature");
        }
    };
    compare_weights(clipped, "absorbed");
    compare_weights(reflected, "reflected");

    // Non-periodic wall stencils must stay inside the stored ghost layer,
    // including a point exactly on each high/low wall.
    Field3d field(d.size(), 3);
    field.data().setOnes();
    for (const Vector3R position : {Vector3R(0, 2, 2), Vector3R(4, 2, 2), Vector3R(2, 0, 2), Vector3R(2, 4, 2),
                                    Vector3R(2, 2, 0), Vector3R(2, 2, 4)}) {
        check::vector(gather_electric(field, position, d), Vector3R(1), 2e-14,
                      "non-periodic wall E stencil stays in range");
        check::vector(gather_magnetic(field, position, d), Vector3R(1), 2e-14,
                      "non-periodic wall B stencil stays in range");
    }
}

void point_wall_simulation() {
    const auto walls = nlohmann::json::parse(R"([
        {"open": [{"face": "XMIN"}, {"face": "XMAX"},
                  {"face": "YMIN"}, {"face": "YMAX"},
                  {"face": "ZMIN"}, {"face": "ZMAX"}]}
    ])");
    for (const std::string assembly : {"reference", "cell"}) {
        auto cfg = settings("ecsim_whitney");
        cfg["WhitneyAssembly"] = assembly;
        cfg["Boundary_conditions"] = walls;
        char program[] = "ecsim_whitney";
        char* argv[] = {program};
        auto sim = build_simulation(cfg, {{"particles", nlohmann::json::array()}}, 1, argv);
        sim->init();
        auto& s = dynamic_cast<SimulationEcsimWhitney&>(*sim);
        auto make = [&](const std::string& name, double charge, double mass) {
            return s.make_particles_array(
                {{"Name", name}, {"Charge", charge}, {"Mass", mass}, {"Density", 1.0}, {"NumPartPerCell", 1}});
        };
        auto electrons = make("Electrons", -1.0, 1.0);
        auto ions = make("Ions", 1.0, 25.0);
        auto neutrals = make("Neutrals", 0.0, 1.0);
        electrons->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
        electrons->add_particle(Particle(Vector3R(2.0, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
        ions->add_particle(Particle(Vector3R(2.0, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
        neutrals->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
        s.species.emplace("Electrons", std::move(electrons));
        s.species.emplace("Ions", std::move(ions));
        s.species.emplace("Neutrals", std::move(neutrals));
        for (auto& [name, sp] : s.species)
            if (!sp->is_neutral())
                s.charged_species.emplace(name, std::ref(*sp));

        s.prepare_step(1);
        s.make_step(1);
        check::require(s.fieldJp.data().allFinite(), "point wall response current must stay finite");
        check::require(std::isfinite(s.mesh.Lmat2.norm()), "point wall response matrix must stay finite");
        check::require(s.species.at("Electrons")->get_total_num_of_particles() == 1,
                       "point scheme must remove only the crossing electron");
        check::require(s.species.at("Neutrals")->get_total_num_of_particles() == 0,
                       "point scheme must remove a neutral crossing the wall");
        check::require(s.species.at("Ions")->get_total_num_of_particles() == 1,
                       "point scheme must keep a particle away from walls");
        s.finalize();
    }
}

void whitney_corr_wall_simulation() {
    const auto walls = nlohmann::json::parse(R"([
        {"open": [{"face": "XMIN"}, {"face": "XMAX"},
                  {"face": "YMIN"}, {"face": "YMAX"},
                  {"face": "ZMIN"}, {"face": "ZMAX"}]}
    ])");
    auto cfg = settings("ecsim_whitney_corr");
    cfg["Boundary_conditions"] = walls;
    char program[] = "ecsim_whitney";
    char* argv[] = {program};
    auto sim = build_simulation(cfg, {{"particles", nlohmann::json::array()}}, 1, argv);
    sim->init();
    auto& s = dynamic_cast<SimulationEcsim&>(*sim);
    auto make = [&](const std::string& name, double charge, double mass) {
        return s.make_particles_array(
            {{"Name", name}, {"Charge", charge}, {"Mass", mass}, {"Density", 1.0}, {"NumPartPerCell", 1}});
    };
    auto electrons = make("Electrons", -1.0, 1.0);
    auto ions = make("Ions", 1.0, 25.0);
    auto neutrals = make("Neutrals", 0.0, 1.0);
    electrons->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
    electrons->add_particle(Particle(Vector3R(2.0, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
    ions->add_particle(Particle(Vector3R(2.0, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
    neutrals->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
    s.species.emplace("Electrons", std::move(electrons));
    s.species.emplace("Ions", std::move(ions));
    s.species.emplace("Neutrals", std::move(neutrals));
    for (auto& [name, sp] : s.species)
        if (!sp->is_neutral())
            s.charged_species.emplace(name, std::ref(*sp));

    s.make_step(1);
    check::require(s.fieldJp.data().allFinite() && s.fieldJe.data().allFinite(),
                   "Whitney-Corr wall currents must stay finite");
    check::require(s.fieldEn.data().allFinite() && s.fieldBn.data().allFinite(),
                   "Whitney-Corr wall fields must stay finite");
    check::require(s.species.at("Electrons")->get_total_num_of_particles() == 1,
                   "Whitney-Corr must remove only the crossing electron");
    check::require(s.species.at("Neutrals")->get_total_num_of_particles() == 0,
                   "Whitney-Corr must remove a neutral crossing the wall");
    check::require(s.species.at("Ions")->get_total_num_of_particles() == 1,
                   "Whitney-Corr must keep a particle away from walls");
    s.finalize();
}

void frozen_ballistic_wall_simulation() {
    const auto walls = nlohmann::json::parse(R"([
        {"open": [{"face": "XMIN"}, {"face": "XMAX"},
                  {"face": "YMIN"}, {"face": "YMAX"},
                  {"face": "ZMIN"}, {"face": "ZMAX"}]}
    ])");
    for (const std::string assembly : {"reference", "cell"}) {
        auto cfg = settings("ecsim_whitney_frozen");
        cfg["WhitneyOrbit"] = {{"predictor", "ballistic"}, {"assembly", assembly}};
        cfg["Boundary_conditions"] = walls;
        char program[] = "ecsim_whitney";
        char* argv[] = {program};
        auto sim = build_simulation(cfg, {{"particles", nlohmann::json::array()}}, 1, argv);
        sim->init();
        auto& s = dynamic_cast<SimulationEcsimWhitneyFrozen&>(*sim);
        auto make = [&](const std::string& name, double charge, double mass) {
            return s.make_particles_array(
                {{"Name", name}, {"Charge", charge}, {"Mass", mass}, {"Density", 1.0}, {"NumPartPerCell", 1}});
        };
        auto electrons = make("Electrons", -1.0, 1.0);
        auto ions = make("Ions", 1.0, 25.0);
        auto neutrals = make("Neutrals", 0.0, 1.0);
        electrons->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
        ions->add_particle(Particle(Vector3R(2.0, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
        neutrals->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
        s.species.emplace("Electrons", std::move(electrons));
        s.species.emplace("Ions", std::move(ions));
        s.species.emplace("Neutrals", std::move(neutrals));
        for (auto& [name, sp] : s.species)
            if (!sp->is_neutral())
                s.charged_species.emplace(name, std::ref(*sp));

        s.prepare_step(1);
        s.make_step(1);
        check::require(s.fieldJp.data().allFinite(), "wall response current must stay finite");
        check::require(std::isfinite(s.mesh.Lmat2.norm()), "wall response matrix must stay finite");
        check::require(s.species.at("Electrons")->get_total_num_of_particles() == 0,
                       "an electron crossing an absorbing wall must be removed");
        check::require(s.species.at("Neutrals")->get_total_num_of_particles() == 0,
                       "a neutral crossing an absorbing wall must be removed");
        check::require(s.species.at("Ions")->get_total_num_of_particles() == 1,
                       "a particle away from walls must survive");
        s.finalize();
    }
}

std::unique_ptr<Simulation> wall_calculation(const nlohmann::json& boundaries, const std::string& assembly) {
    auto cfg = settings("ecsim_whitney_frozen");
    cfg["WhitneyOrbit"] = {{"predictor", "ballistic"}, {"assembly", assembly}};
    cfg["Boundary_conditions"] = boundaries;
    char program[] = "ecsim_whitney";
    char* argv[] = {program};
    auto sim = build_simulation(cfg, {{"particles", nlohmann::json::array()}}, 1, argv);
    sim->init();
    return sim;
}

void add_wall_plasma(SimulationEcsimWhitneyFrozen& s) {
    auto make = [&](const std::string& name, double charge, double mass) {
        return s.make_particles_array(
            {{"Name", name}, {"Charge", charge}, {"Mass", mass}, {"Density", 1.0}, {"NumPartPerCell", 1}});
    };
    auto electrons = make("Electrons", -1.0, 1.0);
    auto ions = make("Ions", 1.0, 25.0);
    auto neutrals = make("Neutrals", 0.0, 1.0);
    // Two crossings exercise both open axes; the second is periodic in the
    // mixed configuration and must wrap instead of being absorbed.
    electrons->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
    electrons->add_particle(Particle(Vector3R(2.0, 2.0, 5.1), Vector3R(0.0, 0.0, 10.0)));
    ions->add_particle(Particle(Vector3R(2.0, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
    ions->add_particle(Particle(Vector3R(0.5, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
    neutrals->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
    s.species.emplace("Electrons", std::move(electrons));
    s.species.emplace("Ions", std::move(ions));
    s.species.emplace("Neutrals", std::move(neutrals));
    for (auto& [name, sp] : s.species)
        if (!sp->is_neutral())
            s.charged_species.emplace(name, std::ref(*sp));
}

void frozen_ballistic_wall_matrix() {
    const auto all_walls = nlohmann::json::parse(R"([
        {"open": [{"face": "XMIN"}, {"face": "XMAX"},
                  {"face": "YMIN"}, {"face": "YMAX"},
                  {"face": "ZMIN"}, {"face": "ZMAX"}]}
    ])");
    // Mirror-trap-like mixed box: periodic X/Y, absorbing Z.
    const auto periodic_xy = nlohmann::json::parse(R"([
        {"periodic": {"face": "XMIN"}},
        {"periodic": {"face": "YMIN"}},
        {"open": [{"face": "ZMIN"}, {"face": "ZMAX"}]}
    ])");
    for (const auto& boundaries : {all_walls, periodic_xy}) {
        auto reference = wall_calculation(boundaries, "reference");
        auto optimized = wall_calculation(boundaries, "cell");
        auto& a = dynamic_cast<SimulationEcsimWhitneyFrozen&>(*reference);
        auto& b = dynamic_cast<SimulationEcsimWhitneyFrozen&>(*optimized);
        add_wall_plasma(a);
        add_wall_plasma(b);
        for (auto* s : {&a, &b}) {
            s->prepare_step(1);
            s->drift_particles_and_assemble_response();
        }
        check::near((a.mesh.Lmat2 - b.mesh.Lmat2).norm(), 0, 2e-14,
                    "wall cell assembly preserves the complete orbit L matrix");
        check::near((a.fieldJp - b.fieldJp).norm(), 0, 2e-14,
                    "wall cell assembly preserves the clipped orbit predictor current");
        for (auto* s : {&a, &b}) {
            for (int i = 0; i < s->fieldEp.capacity(); ++i) {
                s->fieldEp[i] = 0.2 * std::cos(0.11 * i);
                s->fieldE_external[i] = 0.05 * std::sin(0.13 * i);
            }
            ecsimWhitney::sync_periodic(s->fieldEp, s->domain);
            ecsimWhitney::sync_periodic(s->fieldE_external, s->domain);
            s->advance_particle_velocities();
        }
        for (const auto& [name, sp] : a.species) {
            const auto& other = *b.species.at(name);
            check::near((sp->currentOnGrid - other.currentOnGrid).norm(), 0, 2e-13,
                        "wall push preserves the complete adjoint current");
            check::require(sp->get_total_num_of_particles() == other.get_total_num_of_particles(),
                           "wall reference and cell removal agree");
            for (int cell = 0; cell < sp->size(); ++cell) {
                const auto& particles = sp->particlesData(cell);
                const auto& expected = other.particlesData(cell);
                check::require(particles.size() == expected.size(), "wall push preserves particle buckets");
                for (size_t p = 0; p < particles.size(); ++p) {
                    check::vector(particles[p].coord, expected[p].coord, 2e-14,
                                  "wall reference and cell coordinates agree");
                    check::vector(particles[p].velocity, expected[p].velocity, 2e-14,
                                  "wall reference and cell velocities agree");
                }
            }
        }
        a.finalize();
        b.finalize();
    }
}

void frozen_ballistic_wall_energy() {
    const auto walls = nlohmann::json::parse(R"([
        {"open": [{"face": "XMIN"}, {"face": "XMAX"},
                  {"face": "YMIN"}, {"face": "YMAX"},
                  {"face": "ZMIN"}, {"face": "ZMAX"}]}
    ])");
    auto sim = wall_calculation(walls, "reference");
    auto& s = dynamic_cast<SimulationEcsimWhitneyFrozen&>(*sim);
    auto make = [&](const std::string& name, double charge, double mass) {
        return s.make_particles_array(
            {{"Name", name}, {"Charge", charge}, {"Mass", mass}, {"Density", 1.0}, {"NumPartPerCell", 1}});
    };
    auto electrons = make("Electrons", -1.0, 1.0);
    auto ions = make("Ions", 1.0, 25.0);
    electrons->add_particle(Particle(Vector3R(0.2, 2.0, 2.0), Vector3R(-10.0, 0.0, 0.0)));
    ions->add_particle(Particle(Vector3R(2.0, 2.0, 2.0), Vector3R(0.0, 0.0, 0.0)));
    s.species.emplace("Electrons", std::move(electrons));
    s.species.emplace("Ions", std::move(ions));
    for (auto& [name, sp] : s.species)
        if (!sp->is_neutral())
            s.charged_species.emplace(name, std::ref(*sp));
    s.fieldE_external(1, 1, 1, 0) = 0.5;
    s.prepare_step(1);
    s.make_step(1);
    Diagnostics diagnostic(s.species);
    s.diagnostic_energy(diagnostic);
    check::require(s.species.at("Electrons")->get_total_num_of_particles() == 0,
                   "absorbed electron must leave the species");
    check::require(s.species.at("Electrons")->diag.total_lost_energy() > 0,
                   "absorbed electron must contribute loss energy");
    check::near(diagnostic.energy.at("energyConserve"), 0, 1e-11,
                "wall absorption energy balance includes work on the removed orbit");
    s.finalize();
}

void rejects_walls() {
    auto cfg = settings("ecsim_whitney_frozen");
    cfg["WhitneyOrbit"] = {{"predictor", "old_fields"}};
    cfg["Boundary_conditions"] = nlohmann::json::parse(R"([
        {"open": [{"face": "ZMIN"}, {"face": "ZMAX"}]},
        {"periodic": [{"face": "XMIN"}, {"face": "XMAX"},
                      {"face": "YMIN"}, {"face": "YMAX"}]}
    ])");
    char program[] = "ecsim_whitney";
    char* argv[] = {program};
    auto sim = build_simulation(cfg, {{"particles", nlohmann::json::array()}}, 1, argv);
    bool rejected = false;
    try {
        sim->init();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "unsupported walls must fail before the calculation starts");
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("Whitney tensor reference, constants and periodic seams", kernel_properties);
        check::scenario("Compatible CIC scalar gradient and Whitney edge gather", compatible_gradient);
        check::scenario("Two-leg Whitney current preserves CIC charge across cell and periodic planes",
                        two_leg_orbit_continuity);
        check::scenario("Integrated Whitney current preserves charge across cell and periodic planes",
                        orbit_continuity);
        check::scenario("Analytic orbit weights: quadrature, partition, crossings, reversal and invalid paths",
                        analytic_orbit_properties);
        check::scenario("Frozen orbit quadrature, zero-length limit, reversal and constant-field predictors",
                        frozen_orbit_limits_and_predictors);
        check::scenario("Assembled response, midpoint Boris current and electric work",
                        [] { response_matches_boris_work(); });
        check::scenario("Parallel point response, midpoint Boris current and electric work",
                        [] { response_matches_boris_work("", "cell"); });
        check::scenario("Parallel point complete matrix and particle push",
                        [] { cell_assembly_matches_reference(""); });
        for (const std::string predictor : {"ballistic", "old_fields", "half_drift", "centered_fields"}) {
            check::scenario("Frozen " + predictor + ": response, Boris current and electric work",
                            [&] { response_matches_boris_work(predictor); });
            check::scenario("Frozen cell " + predictor + ": response, Boris current and electric work",
                            [&] { response_matches_boris_work(predictor, "cell"); });
            check::scenario("Frozen cell " + predictor + ": complete reference matrix",
                            [&] { cell_assembly_matches_reference(predictor); });
        }
        for (const std::string scheme : {"ecsim_whitney_corr", "ecsim_whitney_frozen_corr"}) {
            check::scenario(scheme + " path current preserves CIC charge in a simulation step",
                            [&] { corr_charge_continuity(scheme); });
            check::scenario(scheme + " conserves total energy over 20 periodic steps",
                            [&] { corr_energy_conservation(scheme); });
        }
        check::scenario("Periodic plasma: energy conservation and ECSIM comparison", compare_schemes);
        check::scenario("Parallel frozen predictor propagates non-finite orbit failure",
                        cell_predictor_failure_propagates);
        check::scenario("Cached orbit response: invalidation and rebuild after injection", cached_response_lifetime);
        check::scenario("Ballistic wall path clipping and non-periodic stencils", ballistic_wall_path_clipping);
        check::scenario("Point wall simulation removes crossing particles", point_wall_simulation);
        check::scenario("Whitney-Corr absorbing wall removes crossing particles", whitney_corr_wall_simulation);
        check::scenario("Frozen ballistic wall simulation absorbs crossing particles",
                        frozen_ballistic_wall_simulation);
        check::scenario("Frozen ballistic wall cell assembly matches reference", frozen_ballistic_wall_matrix);
        check::scenario("Frozen ballistic wall absorption preserves energy balance", frozen_ballistic_wall_energy);
        check::scenario("Unsupported wall configuration is rejected", rejects_walls);
    });
}
