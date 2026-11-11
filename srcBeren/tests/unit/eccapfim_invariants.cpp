#include <filesystem>
#include <fstream>
#include <limits>
#include <random>

#include "simulation/eccapfim/algorithms_eccapfim.h"
#include "simulation/eccapfim/simulation_eccapfim.h"
#include "config/boundary_config.h"
#include "io/field_io.h"
#include "test_support.h"

namespace {
constexpr int Z = static_cast<int>(Axis::Z);
double physical_sum(const Field3d& field, const Domain& domain) {
    const int g = GHOST_CELLS;
    const auto n = domain.num_cells();
    double result = 0;
    for (int i = g; i < g + n.x(); ++i)
        for (int j = g; j < g + n.y(); ++j)
            for (int k = g; k < g + n.z(); ++k)
                for (int c = 0; c < field.nd(); ++c) result += field(i, j, k, c);
    return result;
}

void continuity(const Field3d& before, const Field3d& after, const Field3d& current, double dt, const Domain& domain,
                double tolerance) {
    const int g = GHOST_CELLS;
    const auto n = domain.num_cells();
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k) {
                const Vector3I node(i + g, j + g, k + g);
                double divergence = 0;
                for (int d = 0; d < 3; ++d) {
                    auto prev = node;
                    prev[d] = algorithmsECCAPFIM::periodic_node(node[d] - g - 1, d, domain);
                    divergence +=
                        (current(node.x(), node.y(), node.z(), d) - current(prev.x(), prev.y(), prev.z(), d)) /
                        domain.cell_size(d);
                }
                check::near(
                    after(node.x(), node.y(), node.z(), 0) - before(node.x(), node.y(), node.z(), 0) + dt * divergence,
                    0, tolerance, "pointwise quadratic charge continuity");
            }
}

void orbit_invariants() {
    for (const int ghost : {GHOST_CELLS}) {
        Domain domain;
        domain.init({5, 6, 7}, {0.7, 1.1, 0.4});
        Field3d electric(domain.size(), 3), magnetic(domain.size(), 3);
        magnetic.setZero();
        for (int p = 0; p < electric.capacity(); ++p) electric[p] = std::sin(0.13 * p);
        algorithmsECCAPFIM::copy_periodic_nodes(electric, domain);
        // Zero orbit, negative coordinates, exact half-grid planes, reverse paths,
        // corner crossing and several whole periods. Different x/y/z strides matter.
        const std::vector<std::pair<Vector3R, Vector3R>> paths = {
            {{0.2, 0.3, 0.4}, {0.2, 0.3, 0.4}},       {{-0.35, 0.55, 0.2}, {0.35, -0.55, -0.2}},
            {{3.49, 6.59, 2.79}, {3.71, 6.72, 2.93}}, {{0.8, 1.2, 0.7}, {8.4, -7.5, 5.7}},
            {{8.4, -7.5, 5.7}, {0.8, 1.2, 0.7}},      {{0.7, 1.1, 0.4}, {0.7, 1.1, 1.0}}};
        for (const auto& [start, end] : paths) {
            const double dt = 0.43, charge = -0.7;
            Field3d before(domain.size(), 1), after(domain.size(), 1), current(domain.size(), 3);
            before.setZero();
            after.setZero();
            current.setZero();
            algorithmsECCAPFIM::deposit_charge(before, start, charge, domain);
            algorithmsECCAPFIM::deposit_charge(after, end, charge, domain);
            const auto velocity = (end - start) / dt;
            algorithmsECCAPFIM::deposit_orbit(current, start, end, velocity, charge, domain);
            continuity(before, after, current, dt, domain, 2e-14);
            check::near(physical_sum(before, domain), charge, 2e-15, "charge normalization");
            const auto gather = algorithmsECCAPFIM::gather_orbit(electric, magnetic, start, end, domain);
            // Explicit physical range avoids counting the duplicated periodic nodes.
            IndexRange range;
            range.start = Vector3I(ghost);
            range.end = domain.num_cells() + Vector3I(ghost);
            check::near(dot_product_sum(electric, current, range), charge * velocity.dot(gather.electric), 3e-14,
                        "electric gather/current scatter work identity");
        }
        electric.data().setConstant(0.3);
        magnetic.data().setConstant(-0.2);
        const auto gathered =
            algorithmsECCAPFIM::gather_orbit(electric, magnetic, {-0.5, 0.3, 2.7}, {4.5, 7.0, -0.3}, domain);
        check::vector(gathered.electric, Vector3R(0.3), 2e-15, "constant electric field partition of unity");
        check::vector(gathered.magnetic, Vector3R(-0.2), 2e-15, "quadratic magnetic gather partition of unity");

        // Uniform B gives the analytic Crank-Nicolson Cayley rotation, |v| conserved.
        electric.setZero();
        magnetic.setZero();
        const int g = GHOST_CELLS;
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k) magnetic(i, j, k, 2) = 1.7;
        Field3d current(domain.size(), 3);
        current.setZero();
        const Particle initial({3.49, 6.59, 2.79}, {0.4, -0.2, 0.1});
        const double dt = 2.3, qm = -0.8;
        const auto end = algorithmsECCAPFIM::advance_particle(initial, electric, magnetic, current, qm, -0.7, dt,
                                                              domain, 1e-13, 100);
        const Vector3R b(0, 0, 0.5 * dt * qm * 1.7);
        const auto vh = (initial.velocity + initial.velocity.cross(b) + b * initial.velocity.dot(b)) / (1 + b.dot(b));
        check::vector(end.velocity, 2.0 * vh - initial.velocity, 2e-13, "analytic CN gyro rotation");
        check::near(end.velocity.dot(end.velocity), initial.velocity.dot(initial.velocity), 2e-14,
                    "magnetic work is zero");
        check::vector(end.coord, algorithmsECCAPFIM::wrap_position(initial.coord + dt * vh, domain), 2e-13,
                      "CN midpoint drift");
        check::require(algorithmsECCAPFIM::periodic_node(-1, 0, domain) == g + 4, "negative physical index mapping");

        // A spatially varying force with only four fixed-point iterations forces
        // time subdivision. Its J must still use the full dt normalization.
        for (int p = 0; p < electric.capacity(); ++p) electric[p] = 0.2 * std::sin(0.13 * p);
        algorithmsECCAPFIM::copy_periodic_nodes(electric, domain);
        magnetic.setZero();
        current.setZero();
        const auto subcycled =
            algorithmsECCAPFIM::advance_particle(initial, electric, magnetic, current, qm, -0.7, dt, domain, 1e-10, 4);
        Field3d before(domain.size(), 1), after(domain.size(), 1);
        before.setZero();
        after.setZero();
        algorithmsECCAPFIM::deposit_charge(before, initial.coord, -0.7, domain);
        algorithmsECCAPFIM::deposit_charge(after, subcycled.coord, -0.7, domain);
        continuity(before, after, current, dt, domain, 3e-13);
        IndexRange range;
        range.start = Vector3I(ghost);
        range.end = domain.num_cells() + Vector3I(ghost);
        const double kinetic_change = 0.5 * (-0.7 / qm) * (subcycled.velocity.squaredNorm() - initial.velocity.squaredNorm());
        check::near(kinetic_change, dt * dot_product_sum(electric, current, range), 2e-10,
                    "substep particle energy equals full-step current work");
    }
}

nlohmann::json settings() {
    return {{"Scheme", "eccapfim"},
            {"Dt", 0.35},
            {"NumCellsX", 4},
            {"NumCellsY", 5},
            {"NumCellsZ", 3},
            {"Dx", 0.7},
            {"Dy", 1.1},
            {"Dz", 0.9},
            {"StartFromTime", 0},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

std::unique_ptr<SimulationEccapfim> simulation(nlohmann::json config, bool particles = true) {
    auto cfg = nlohmann::json{{"particles", nlohmann::json::array()}};
    if (particles)
        cfg["particles"] =
            nlohmann::json::array({check::species("Electrons", -1, 1, 1, 1), check::species("Ions", 1, 9, 1, 1),
                                   check::species("Gas", 0, 2, 1, 1)});
    char name[] = "eccapfim_invariants";
    char* argv[] = {name};
    auto base = build_simulation(config, cfg, 1, argv);
    auto* concrete = dynamic_cast<SimulationEccapfim*>(base.get());
    check::require(concrete != nullptr, "scheme factory dispatch");
    base.release();
    std::unique_ptr<SimulationEccapfim> result(concrete);
    result->init();
    return result;
}

double field_energy(const Field3d& field, const IndexRange& range) {
    return 0.5 * dot_product_sum(field, field, range);
}

double energy(const SimulationEccapfim& sim) {
    const auto range = sim.bc_handler.active_range(sim.domain.grid);
    double result = field_energy(sim.fieldEn, range) + field_energy(sim.fieldBn, range);
    for (const auto& [name, sp] : sim.species) result += sp->get_kinetic_energy();
    return result;
}

nlohmann::json wall_settings(bool reflecting = false) {
    auto config = settings();
    config["Boundary_conditions"] = nlohmann::json::parse(
        R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},
             {"open":[{"face":"ZMIN"},{"face":"ZMAX"}]}])");
    if (reflecting)
        config["Boundary_conditions"][2] = {
            {"electron_reflection",
             nlohmann::json::array({{{"face", "ZMIN"}, {"radius", 100}, {"energy_threshold", 1e6}},
                                    {{"face", "ZMAX"}, {"radius", 100}, {"energy_threshold", 1e6}}})}};
    return config;
}

void wall_orbits() {
    for (bool reflecting : {false, true}) {
        auto sim = simulation(wall_settings(reflecting));
        auto& sp = *sim->species.at("Electrons");
        const auto& domain = sim->domain;
        const auto policy = sim->bc_handler.path_policy(sp, domain);
        Field3d electric(domain.size(), 3), magnetic(domain.size(), 3), current(domain.size(), 3);
        electric.setZero();
        magnetic.setZero();
        const double length = domain.geom.box_max.z(), dt = 0.8;
        for (bool high : {false, true}) {
            current.setZero();
            const double sign = high ? 1 : -1;
            const Particle initial({2.79, 5.49, high ? length - 0.1 : 0.1}, {0.3, 0.2, sign * 0.5});
            algorithmsECCAPFIM::OrbitBoundaryEvents event;
            const auto end = algorithmsECCAPFIM::advance_particle(initial, electric, magnetic, current, -1, -0.7, dt,
                                                                  domain, 1e-13, 100, &sim->bc_handler, policy, &event);
            const auto face = high ? Face::ZMAX : Face::ZMIN;
            check::near(end.velocity.squaredNorm(), initial.velocity.squaredNorm(), 1e-14, "wall does no kinetic work");
            if (reflecting) {
                check::require(event.removed == Face::UNAVAILABLE, "reflected electron removed");
                check::require(event.reflected[static_cast<size_t>(face)] == 1, "wall reflection count");
                check::near(end.coord.z(), high ? length - 0.3 : 0.3, 1e-13, "remaining-time reflected drift");
                check::near(end.velocity.z(), -initial.velocity.z(), 1e-14, "specular reflected velocity");
            } else {
                check::require(event.removed == face, "absorbed particle fate");
                check::near(end.coord.z(), high ? length : 0, 1e-14, "absorption at exact wall");
                check::near(end.coord.x(), std::fmod(initial.coord.x() + 0.2 * initial.velocity.x(), 2.8), 1e-13,
                            "absorbed orbit ends at hit time, including periodic x seam");
            }
            // The z-integrated current equals q times the actual net displacement / full dt.
            const auto range = sim->bc_handler.active_range(domain.grid);
            double jz = 0;
            for (int i = range.start.x(); i < range.end.x(); ++i)
                for (int j = range.start.y(); j < range.end.y(); ++j)
                    for (int k = range.start.z(); k < range.end.z(); ++k) jz += current(i, j, k, Z);
            check::near(dt * jz, -0.7 * (end.coord.z() - initial.coord.z()), 1e-13,
                        "wall orbit uses full-step current normalization");
        }
        if (reflecting) {
            current.setZero();
            const Particle initial({1.4, 2.75, 0.2}, {0, 0, 1});
            algorithmsECCAPFIM::OrbitBoundaryEvents event;
            const auto end =
                algorithmsECCAPFIM::advance_particle(initial, electric, magnetic, current, -1, -0.7, 3 * length, domain,
                                                     1e-13, 100, &sim->bc_handler, policy, &event);
            check::near(end.coord.z(), length - 0.2, 2e-12, "multiple wall hits in one step");
            check::require(event.reflected[static_cast<size_t>(Face::ZMIN)] == 1 &&
                               event.reflected[static_cast<size_t>(Face::ZMAX)] == 2,
                           "multiple reflection counts");
        }
        sim->finalize();
    }
}

void accelerated_wall_orbits() {
    for (bool reflecting : {false, true}) {
        auto sim = simulation(wall_settings(reflecting));
        const auto& domain = sim->domain;
        const auto policy = sim->bc_handler.path_policy(*sim->species.at("Electrons"), domain);
        Field3d electric(domain.size(), 3), magnetic(domain.size(), 3), current(domain.size(), 3);
        magnetic.setZero();
        for (double speed : {0.0, 0.2}) {
            for (bool high : {false, true}) {
                electric.setZero();
                current.setZero();
                const double sign = high ? 1 : -1;
                const double acceleration = 0.7, distance = 0.1, dt = 0.8;
                for (int i = 0; i < domain.size().x(); ++i)
                    for (int j = 0; j < domain.size().y(); ++j)
                        for (int k = 0; k < domain.size().z(); ++k) electric(i, j, k, Z) = -sign * acceleration;
                const double wall = high ? domain.geom.box_max.z() : 0;
                const Particle initial({1.4, 2.75, wall - sign * distance}, {0, 0, sign * speed});
                const double hit_time = (std::sqrt(speed * speed + 2 * acceleration * distance) - speed) / acceleration;
                const double hit_speed = speed + acceleration * hit_time;
                algorithmsECCAPFIM::OrbitBoundaryEvents event;
                const auto end =
                    algorithmsECCAPFIM::advance_particle(initial, electric, magnetic, current, -1, -0.7, dt, domain,
                                                         1e-12, 100, &sim->bc_handler, policy, &event);
                const double rest = dt - hit_time;
                check::near(end.coord.z(),
                            reflecting ? wall + sign * (-hit_speed * rest + 0.5 * acceleration * rest * rest) : wall,
                            2e-11, "accelerated orbit uses physical wall time");
                check::near(end.velocity.z(), sign * (reflecting ? -hit_speed + acceleration * rest : hit_speed), 2e-11,
                            "velocity at absorption or after reflected remaining time");
                const auto range = sim->bc_handler.active_range(domain.grid);
                check::near(0.5 * 0.7 * (end.velocity.squaredNorm() - initial.velocity.squaredNorm()),
                            dt * dot_product_sum(electric, current, range), 2e-11,
                            "wall orbit kinetic change equals electric current work");
            }
        }
        sim->finalize();
    }
}

void decelerating_wall_orbits() {
    // A weak impact close to a turning point is ill-conditioned for the
    // chord-based iteration tau <- tau * distance / displacement.
    // Uniform opposing E gives an independent analytic time and impact speed.
    const double speed = 1, acceleration = 1, distance = 0.49995, dt = 1;
    const double hit_speed = std::sqrt(speed * speed - 2 * acceleration * distance);
    const double hit_time = 2 * distance / (speed + hit_speed);
    for (bool reflecting : {false, true}) {
        auto sim = simulation(wall_settings(reflecting));
        const auto& domain = sim->domain;
        const auto policy = sim->bc_handler.path_policy(*sim->species.at("Electrons"), domain);
        Field3d electric(domain.size(), 3), magnetic(domain.size(), 3), current(domain.size(), 3);
        magnetic.setZero();
        for (bool high : {false, true}) {
            electric.setZero();
            current.setZero();
            const double sign = high ? 1 : -1;
            const double wall = high ? domain.geom.box_max.z() : domain.geom.box_min.z();
            for (int i = 0; i < domain.size().x(); ++i)
                for (int j = 0; j < domain.size().y(); ++j)
                    for (int k = 0; k < domain.size().z(); ++k) electric(i, j, k, Z) = sign * acceleration;
            const Particle initial({1.4, 2.75, wall - sign * distance}, {0.03, -0.01, sign * speed});
            algorithmsECCAPFIM::OrbitBoundaryEvents event;
            const auto end = algorithmsECCAPFIM::advance_particle(initial, electric, magnetic, current, -1, -0.7, dt,
                                                                  domain, 1e-14, 100, &sim->bc_handler, policy, &event);
            const double rest = dt - hit_time;
            const double final_speed = reflecting ? -hit_speed - acceleration * rest : hit_speed;
            const double final_position =
                reflecting ? wall - sign * (hit_speed * rest + 0.5 * acceleration * rest * rest) : wall;
            const auto face = high ? Face::ZMAX : Face::ZMIN;
            check::near(end.coord.z(), final_position, 1e-12, "decelerating wall endpoint");
            check::near(end.velocity.z(), sign * final_speed, 2e-11, "decelerating wall impact speed");
            check::near(end.coord.x(), initial.coord.x() + initial.velocity.x() * (reflecting ? dt : hit_time), 1e-12,
                        "decelerating wall impact time from transverse drift");
            check::require(event.removed == (reflecting ? Face::UNAVAILABLE : face), "decelerating wall fate");
            check::require(event.reflected[static_cast<size_t>(face)] == (reflecting ? 1 : 0),
                           "decelerating wall reflection count");
            const auto range = sim->bc_handler.active_range(domain.grid);
            check::near(0.5 * 0.7 * (end.velocity.squaredNorm() - initial.velocity.squaredNorm()),
                        dt * dot_product_sum(electric, current, range), 2e-13,
                        "decelerating wall kinetic change equals deposited electric work");
            // Keep the endpoint charge at the wall before applying the
            // absorbing sink; reflected paths retain their endpoint charge.
            Field3d before(domain.size(), 1), after(domain.size(), 1);
            before.setZero();
            after.setZero();
            algorithmsECCAPFIM::deposit_charge(before, initial.coord, -0.7, domain, &sim->bc_handler);
            algorithmsECCAPFIM::deposit_charge(after, end.coord, -0.7, domain, &sim->bc_handler);
            const int ghost = GHOST_CELLS;
            auto nodes = domain.num_cells();
            ++nodes[Z];
            for (int i = 0; i < nodes.x(); ++i)
                for (int j = 0; j < nodes.y(); ++j)
                    for (int k = 0; k < nodes.z(); ++k) {
                        const auto node = Vector3I(i, j, k) + Vector3I(ghost);
                        double divergence = 0;
                        for (int d = 0; d < 3; ++d) {
                            auto prev = node;
                            prev[d] = algorithmsECCAPFIM::node_owner(node[d] - ghost - 1, d, domain, &sim->bc_handler);
                            divergence +=
                                (current(node.x(), node.y(), node.z(), d) - current(prev.x(), prev.y(), prev.z(), d)) /
                                domain.cell_size(d);
                        }
                        check::near(after(node.x(), node.y(), node.z(), 0) - before(node.x(), node.y(), node.z(), 0) +
                                        dt * divergence,
                                    0, 2e-13, "decelerating wall pointwise charge continuity");
                    }
        }
        sim->finalize();
    }
}

void mixed_boundary_deposit() {
    for (const int ghost : {GHOST_CELLS})
        for (int axis = 0; axis < 3; ++axis) {
            Domain domain;
            domain.init({4, 5, 6}, {0.7, 1.1, 0.9});
            const char* low[] = {"XMIN", "YMIN", "ZMIN"};
            const char* high[] = {"XMAX", "YMAX", "ZMAX"};
            auto config = settings();
            config["Boundary_conditions"] = nlohmann::json::array();
            for (int d = 0; d < 3; ++d)
                if (d != axis)
                    config["Boundary_conditions"].push_back({{"periodic", {{"face", low[d]}}}});
            config["Boundary_conditions"].push_back(
                {{"open", nlohmann::json::array({{{"face", low[axis]}}, {{"face", high[axis]}}})}});
            BoundaryConditionHandler bc = boundary_conditions_from_config(config, domain);
            Field3d electric(domain.size(), 3), magnetic(domain.size(), 3);
            magnetic.setZero();
            for (int p = 0; p < electric.capacity(); ++p) electric[p] = std::sin(0.13 * p);
            algorithmsECCAPFIM::copy_periodic_nodes(electric, domain, &bc);
            for (bool upper : {false, true}) {
                Vector3R start = domain.geom.box_max - Vector3R(0.01);
                start[axis] = upper ? domain.geom.box_max[axis] - 0.1 : 0.1;
                Vector3R end = start + Vector3R(0.2);
                end[axis] = upper ? domain.geom.box_max[axis] : 0;
                const double dt = 0.8, charge = -0.7;
                Field3d before(domain.size(), 1), at_wall(domain.size(), 1), current(domain.size(), 3);
                before.setZero();
                at_wall.setZero();
                current.setZero();
                algorithmsECCAPFIM::deposit_charge(before, start, charge, domain, &bc);
                algorithmsECCAPFIM::deposit_charge(at_wall, end, charge, domain, &bc);
                algorithmsECCAPFIM::deposit_orbit(current, start, end, (end - start) / dt, charge, domain, &bc);
                const auto fields = algorithmsECCAPFIM::gather_orbit(electric, magnetic, start, end, domain, &bc);
                check::near(dt * dot_product_sum(electric, current, bc.active_range(domain.grid)),
                            charge * (end - start).dot(fields.electric), 2e-13,
                            "mixed boundary gather/scatter work with periodic corner crossing");
                auto n = domain.num_cells();
                ++n[axis];   // Include both wall charge nodes on the nonperiodic axis.
                for (int i = 0; i < n.x(); ++i)
                    for (int j = 0; j < n.y(); ++j)
                        for (int k = 0; k < n.z(); ++k) {
                            const Vector3I node = Vector3I(i, j, k) + Vector3I(ghost);
                            double divergence = 0;
                            for (int d = 0; d < 3; ++d) {
                                auto prev = node;
                                prev[d] = algorithmsECCAPFIM::node_owner(node[d] - ghost - 1, d, domain, &bc);
                                divergence += (current(node.x(), node.y(), node.z(), d) -
                                               current(prev.x(), prev.y(), prev.z(), d)) /
                                              domain.cell_size(d);
                            }
                            // Before deletion, charge transported to the wall obeys continuity.
                            // Deletion then removes this endpoint shape as a boundary charge sink.
                            check::near(at_wall(node.x(), node.y(), node.z(), 0) -
                                            before(node.x(), node.y(), node.z(), 0) + dt * divergence,
                                        0, 2e-13, "mixed boundary pointwise charge continuity to wall");
                        }
            }
        }
}

void wall_coupled_invariants() {
    for (bool reflecting : {false, true}) {
        auto config = wall_settings(reflecting);
        config["Dt"] = 0.8;
        auto sim = simulation(config);
        // Include a charged and a neutral loss; the reflected electron must
        // remain while ions and neutrals are absorbed by the same classifier.
        sim->species.at("Electrons")->add_particle(Particle({1.4, 2.75, 2.6}, {0.03, 0, 0.4}));
        sim->species.at("Ions")->add_particle(Particle({1.4, 2.75, 2.6}, {-0.01, 0, 0.4}));
        sim->species.at("Gas")->add_particle(Particle({1.4, 2.75, 0.1}, {0, 0, -0.5}));
        const double before = energy(*sim);
        sim->prepare_step(1);
        sim->make_step(1);
        double lost = 0;
        for (const auto& [name, sp] : sim->species) lost += sp->step_state.total_lost_energy();
        check::near(energy(*sim) - before + lost, 0, 3e-9, "coupled energy including wall losses");
        check::require(sim->species.at("Electrons")->get_total_num_of_particles() == (reflecting ? 1 : 0),
                       "accepted electron wall fate");
        check::require(sim->species.at("Ions")->get_total_num_of_particles() == 0, "ion wall loss");
        check::require(sim->species.at("Gas")->get_total_num_of_particles() == 0, "neutral wall loss");
        check::require(sim->species.at("Ions")->step_state.total_lost_count_z() == 1, "iterations duplicated losses");
        if (reflecting)
            check::require(sim->species.at("Electrons")->step_state.boundary.at(Face::ZMAX).reflected_count == 1,
                           "iterations duplicated reflections");
        // Every constrained component uses its Yee position, including ghost layers.
        auto masked = sim->fieldEn;
        sim->bc_handler.apply_to_fields(masked, FieldType::ELECTRIC, sim->domain);
        check::near((masked - sim->fieldEn).norm(), 0, 0, "accepted field violates wall mask");
        sim->make_diagnostic(1);
        sim->prepare_step(2);
        check::require(sim->species.at("Ions")->step_state.total_lost_count_z() == 0, "loss counters not reset");
        sim->finalize();
    }
}

void coupled_invariants() {
    for (bool particles : {false, true}) {
        auto sim = simulation(settings(), particles);
        const auto size = sim->domain.size();
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k) {
                    sim->fieldEn(i, j, k, 1) = 0.02 * std::sin(2 * M_PI * (i - 1) / 4);
                    sim->fieldBn(i, j, k, 2) = 0.03 * std::cos(2 * M_PI * (i - 1) / 4);
                }
        algorithmsECCAPFIM::copy_periodic_nodes(sim->fieldEn, sim->domain);
        algorithmsECCAPFIM::copy_periodic_nodes(sim->fieldBn, sim->domain);
        if (particles) {
            // Charge-neutral initial positions include a corner and the periodic seam.
            for (const auto& r : {Vector3R(2.79, 5.49, 2.69), Vector3R(0.35, 1.1, 0.45), Vector3R(1.5, 2.3, 1.4)}) {
                sim->species.at("Electrons")->add_particle(Particle(r, {0.7, -0.3, 0.2}));
                sim->species.at("Ions")->add_particle(Particle(r, {-0.1, 0.2, 0.05}));
            }
            sim->species.at("Gas")->add_particle(Particle({2.79, 5.49, 2.69}, {5, -4, 3}));
        }
        const double initial = energy(*sim);
        for (int step = 1; step <= 5; ++step) {
            Field3d before(sim->domain.size(), 1), after(sim->domain.size(), 1);
            sim->collect_charge_density(before);
            sim->prepare_step(step);
            sim->make_step(step);
            sim->collect_charge_density(after);
            continuity(before, after, sim->fieldJ, sim->dt_, sim->domain, 3e-12);
            Field3d gauss(sim->domain.size(), 1);
            spmv(sim->maxwell.divE, sim->fieldEn, gauss);
            gauss.data() -= after.data();
            check::near(gauss.norm(), 0, 3e-9, "Gauss law propagated from initially neutral plasma");
            check::near(energy(*sim), initial, 2e-9, "coupled particle+Maxwell energy conservation");
            check::require(std::isfinite(sim->nonlinear_residual) && sim->nonlinear_residual < 1e-9,
                           "accepted field residual");
            sim->make_diagnostic(step);
        }
        if (particles) {
            const auto& gas = *sim->species.at("Gas");
            check::require(gas.get_total_num_of_particles() == 1, "long neutral orbit rebin lost a particle");
            for (int c = 0; c < gas.size(); ++c)
                for (const auto& p : gas(c)) check::vector(p.velocity, {5, -4, 3}, 0, "neutral accelerated by fields");
        }
        sim->finalize();
        sim->finalize();
    }
}

void cold_plasma() {
    auto config = settings();
    config["NumCellsZ"] = 1;   // Degenerate periodic direction, including all duplicated stencil columns.
    config["Dt"] = 3.0;        // Resolves an implicit oscillation beyond the explicit plasma frequency limit.
    auto sim = simulation(config);
    const int g = GHOST_CELLS;
    const auto n = sim->domain.num_cells();
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k) {
                const Vector3R r((i + 0.5) * 0.7, (j + 0.5) * 1.1, (k + 0.5) * 0.9);
                sim->species.at("Electrons")->add_particle(Particle(r, {0.01, 0, 0}));
                sim->species.at("Ions")->add_particle(Particle(r, {0, 0, 0}));
            }
    const double initial = energy(*sim);
    sim->prepare_step(1);
    sim->make_step(1);
    const double dt = sim->dt_, omega2 = 1 + 1.0 / 9;
    const double eh = 0.5 * dt * 0.01 / (1 + 0.25 * dt * dt * omega2);
    check::near(sim->fieldEp(g, g, g, 0), eh, 2e-11, "uniform plasma analytic implicit field");
    for (int cell = 0; cell < sim->species.at("Electrons")->size(); ++cell)
        for (const auto& p : (*sim->species.at("Electrons"))(cell))
            check::near(p.velocity.x(), 0.01 - dt * eh, 7e-11, "uniform plasma analytic electron velocity");
    check::near(energy(*sim), initial, 2e-9, "large dt plasma energy conservation");
    sim->finalize();
}

// The mass matrix is (dt/2) dJ/dE of frozen orbits. For a particle at rest
// the orbit does not move, so each column equals the finite-difference
// current response to a unit E at one Yee edge.
struct MassProbe : SimulationEccapfim {
    using SimulationEccapfim::SimulationEccapfim;
    Vector current(const Field3d& electric) {
        evaluate_current(electric + fieldE_external, fieldB + fieldBInit);
        return pack(fieldJ);
    }
    void start() {
        kinetic_states_ = save_particles();
    }
    Eigen::SparseMatrix<double> mass() {
        rebuild_mass_matrix(fieldB + fieldBInit);
        return mass_only_matrix();
    }
    // Solves A x = rhs with the given GMRES preconditioner; returns iterations.
    int solve(YeeStencilMatrix::Preconditioner preconditioner, const Vector& rhs, Vector& x) {
        stencil_operator_.factor(preconditioner);
        int iterations = 0;
        check::require(stencil_operator_.solve_gmres(rhs, x, 1e-10, 50, 2000, iterations), "mass-matrix GMRES");
        return iterations;
    }
    void reset_operator() {
        stencil_operator_.reset();
    }

    int rows() const {
        return stencil_operator_.rows();
    }

    bool circulant() const {
        return stencil_operator_.supports_circulant();
    }
    int dof(int i, int j, int k, int c) const {
        return reduced_dof(i + GHOST_CELLS, j + GHOST_CELLS, k + GHOST_CELLS, c);
    }
};

void mass_matrix_response() {
    auto config = settings();
    config["NumCellsX"] = 8;
    config["NumCellsY"] = 8;
    config["NumCellsZ"] = 8;
    config["Dt"] = 1.0;
    config["Eccapfim"] = {{"preconditioner", "mass_matrix"}, {"particle_tolerance", 1e-14}};
    const auto cfg = nlohmann::json{{"particles", nlohmann::json::array({check::species("Electrons", -1, 1, 1, 1)})}};
    char name[] = "eccapfim_invariants";
    char* argv[] = {name};
    MassProbe sim(config, cfg, 1, argv);
    sim.init();
    for (int p = 0; p < sim.fieldBInit.capacity(); p += 3) {
        sim.fieldBInit[p] = 0.3;
        sim.fieldBInit[p + 1] = -0.2;
        sim.fieldBInit[p + 2] = 0.7;
    }
    sim.species.at("Electrons")->add_particle(Particle({3.3, 4.1, 2.7}, {0, 0, 0}));
    sim.start();
    const auto mass = sim.mass();
    Field3d electric(sim.domain.size(), 3);
    electric.setZero();
    const auto base = sim.current(electric);
    const double eps = 1e-6;
    for (int c = 0; c < 3; ++c)
        for (const Vector3I node : {Vector3I(3, 4, 3), Vector3I(2, 4, 3), Vector3I(4, 5, 2)}) {
            electric.setZero();
            electric(node.x() + GHOST_CELLS, node.y() + GHOST_CELLS, node.z() + GHOST_CELLS, c) = eps;
            algorithmsECCAPFIM::copy_periodic_nodes(electric, sim.domain);
            const auto response = ((0.5 * sim.dt_ / eps) * (sim.current(electric) - base)).eval();
            const Eigen::VectorXd column = mass.col(sim.dof(node.x(), node.y(), node.z(), c));
            check::near((response - column).norm(), 0, 1e-8 * std::max(1.0, response.norm()),
                        "mass-matrix column equals the frozen-orbit current response");
        }
    sim.finalize();
}

// `cells`: grid size; one-cell axes (a column along z) have a trivial
// Fourier transform on their axis. `walls`: conducting (open) faces on z,
// handled by the image extension of the circulant preconditioner.
void mass_matrix_preconditioners(const Vector3I& cells, bool walls = false) {
    auto config = settings();
    config["NumCellsX"] = cells.x();
    config["NumCellsY"] = cells.y();
    config["NumCellsZ"] = cells.z();
    if (walls)
        config["Boundary_conditions"] = nlohmann::json::parse(
            R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"open":[{"face":"ZMIN"},{"face":"ZMAX"}]}])");
    config["Dt"] = 3.0;
    config["Eccapfim"] = {{"preconditioner", "mass_matrix"}};
    const auto cfg = nlohmann::json{{"particles", nlohmann::json::array({check::species("Electrons", -1, 1, 1, 1)})}};
    char name[] = "eccapfim_invariants";
    char* argv[] = {name};
    MassProbe sim(config, cfg, 1, argv);
    sim.init();
    using P = YeeStencilMatrix::Preconditioner;
    check::require(sim.circulant(), "circulant preconditioner is available");
    Eigen::VectorXd rhs = Eigen::VectorXd::LinSpaced(sim.rows(), -1, 1).array().sin();
    Eigen::VectorXd x;
    // The vacuum operator I + dt^2/4 curl curl is translation invariant: FFT inverts it exactly.
    sim.reset_operator();
    check::require(sim.solve(P::Circulant, rhs, x) <= 1, "circulant inverse is exact for a uniform operator");
    // A nonuniform plasma: every preconditioner gives the same solution.
    // Markers at fixed fractions of the box.
    const Vector3R box(cells.x() * 0.7, cells.y() * 1.1, cells.z() * 0.9);
    for (const auto& f : {Vector3R(0.07, 0.08, 0.07), Vector3R(0.5, 0.74, 0.7), Vector3R(0.93, 0.24, 0.03)})
        sim.species.at("Electrons")
            ->add_particle(Particle(Vector3R(f.x() * box.x(), f.y() * box.y(), f.z() * box.z()), {0, 0, 0}));
    sim.start();
    sim.mass();
    Eigen::VectorXd reference;
    sim.solve(P::BlockJacobi, rhs, reference);
    for (const auto preconditioner : {P::Circulant, P::CirculantBlockJacobi}) {
        sim.solve(preconditioner, rhs, x);
        check::near((x - reference).norm() / reference.norm(), 0, 1e-8, "preconditioners give one solution");
    }
    sim.finalize();
}

void maxwell_preconditioners() {
    for (const std::string preconditioner : {"incomplete_cholesky", "amgcl", "block_jacobi"}) {
        auto config = settings();
        config["Eccapfim"] = {{"preconditioner", preconditioner}};
        auto sim = simulation(config);
        const Vector3R position(1.0, 2.0, 1.5);
        sim->species.at("Electrons")->add_particle(Particle(position, {0.01, 0, 0}));
        sim->species.at("Ions")->add_particle(Particle(position, {0, 0, 0}));
        sim->prepare_step(1);
        sim->make_step(1);
        check::require(sim->fieldEn.data().allFinite() && sim->fieldBn.data().allFinite(),
                       preconditioner + " Maxwell solve stays finite");
        sim->finalize();
    }
}

void failure_is_not_a_step() {
    auto config = settings();
    config["Eccapfim"] = {{"max_iterations", 1}, {"mixing", 0.1}, {"atol", 1e-16}, {"rtol", 0}};
    auto sim = simulation(config);
    const Particle initial({2.79, 1.0, 1.0}, {0.3, 0.1, 0});
    sim->species.at("Electrons")->add_particle(initial);
    sim->prepare_step(1);
    const Field3d previous_e = sim->fieldEn, previous_b = sim->fieldBn;
    bool failed = false;
    try {
        sim->make_step(1);
    } catch (const std::runtime_error&) {
        failed = true;
    }
    check::require(failed, "unconverged nonlinear solve was accepted");
    check::near((sim->fieldEn - previous_e).norm(), 0, 0, "failed solve committed E");
    check::near((sim->fieldBn - previous_b).norm(), 0, 0, "failed solve committed B");
    for (int cell = 0; cell < sim->species.at("Electrons")->size(); ++cell)
        for (const auto& p : (*sim->species.at("Electrons"))(cell)) {
            check::vector(p.coord, initial.coord, 0, "failed solve committed position");
            check::vector(p.velocity, initial.velocity, 0, "failed solve committed velocity");
        }
    sim->finalize();
    config = settings();
    config["Boundary_conditions"] = nlohmann::json::array();
    failed = false;
    try {
        auto unsupported = simulation(config);
    } catch (const std::invalid_argument&) {
        failed = true;
    }
    check::require(failed, "nonperiodic boundaries silently accepted");
}
void external_work() {
    auto sim = simulation(settings());
    sim->fieldE_external.data().setConstant(0.07);
    sim->fieldBInit.data().setConstant(0.3);
    sim->species.at("Electrons")->add_particle(Particle({0.69, 1.09, 0.89}, {0.3, -0.2, 0.1}));
    const double initial = energy(*sim);
    sim->prepare_step(1);
    sim->make_step(1);
    const double work =
        sim->dt_ * dot_product_sum(sim->fieldE_external, sim->fieldJ, sim->bc_handler.active_range(sim->domain.grid));
    check::require(std::abs(work) > 1e-6, "external work scenario has zero signal");
    check::near(energy(*sim) - initial, work, 2e-10, "external electric work enters the energy balance");
    sim->finalize();
}

void particle_failure_dump() {
    auto config = settings();
    config["Eccapfim"] = {{"failure_dump", true}, {"iteration_log", true}};
    auto sim = simulation(config);
    const Particle initial({0.4, 0.6, 0.8}, {0.1, -0.2, 0.3});
    auto& electrons = *sim->species.at("Electrons");
    electrons.add_particle(initial);
    // Deliberate nonfinite Lorentz force must stop before accepting any particles.
    sim->fieldBInit.data().setConstant(std::numeric_limits<double>::quiet_NaN());
    sim->prepare_step(7);
    bool failed = false;
    try {
        sim->make_step(7);
    } catch (const algorithmsECCAPFIM::ParticleSolveFailure& error) {
        failed = true;
        check::require(error.details.at("reason") == "nonfinite_iterate", "specific particle failure reason");
        check::require(error.details.at("retries").size() == 1, "failure includes retry history");
    }
    check::require(failed, "nonfinite particle solve was accepted");
    const auto directory = std::filesystem::path("EccapfimFailure/step_7");
    std::ifstream file(directory / "failure.json");
    nlohmann::json report;
    file >> report;
    check::require(report.at("step") == 7 && report.at("field_iteration") == 0, "failure time and field iteration");
    check::require(report.at("species") == "Electrons", "failure species");
    const int cell = report.at("cell");
    check::vector(electrons(cell)[0].coord, initial.coord, 0, "failed solve changed particle position");
    check::vector(electrons(cell)[0].velocity, initial.velocity, 0, "failed solve changed particle velocity");
    Field3d electric(sim->domain.size(), 3);
    read_field_from_file((directory / "electric.backup").string(), electric);
    check::near((electric - sim->fieldEp - sim->fieldE_external).norm(), 0, 0, "failure electric field roundtrip");
    check::require(std::filesystem::file_size(directory / "magnetic.backup") ==
                       3 * sizeof(int) + sim->fieldBInit.capacity() * sizeof(double),
                   "full failure magnetic field file");
    sim->finalize();
}

void restart() {
    auto config = settings();
    config["diagnostics"]["outputs"] = nlohmann::json::array({{{"type", "recovery"}, {"interval", 1}}});
    auto original = simulation(config);
    original->species.at("Electrons")->add_particle(Particle({2.79, 1.09, 0.89}, {0.3, -0.2, 0.1}));
    original->prepare_step(1);
    original->make_step(1);
    original->make_diagnostic(1);
    check::require(std::filesystem::exists("Recovery/Fields/FieldE.backup"), "actual recovery field file missing");
    config["StartFromTime"] = 1;
    config["diagnostics"]["outputs"] = nlohmann::json::array();
    const auto directory = std::filesystem::current_path();
    std::filesystem::create_directories("restart");
    std::filesystem::current_path("restart");   // The existing reader uses ../Recovery.
    auto resumed = simulation(config);
    std::filesystem::current_path(directory);
    check::near((resumed->fieldEn - original->fieldEn).norm(), 0, 0, "restart E integer time level");
    original->prepare_step(2);
    original->make_step(2);
    resumed->prepare_step(2);
    resumed->make_step(2);
    check::near((resumed->fieldEn - original->fieldEn).norm(), 0, 2e-13, "restart next electric field");
    check::near((resumed->fieldBn - original->fieldBn).norm(), 0, 2e-13, "restart next magnetic field");
    check::near(energy(*resumed), energy(*original), 2e-13, "restart total energy");
    const auto& expected = *original->species.at("Electrons");
    const auto& actual = *resumed->species.at("Electrons");
    for (int cell = 0; cell < expected.size(); ++cell) {
        check::require(actual(cell).size() == expected(cell).size(), "restart particle count");
        for (size_t p = 0; p < expected(cell).size(); ++p) {
            check::vector(actual(cell)[p].coord, expected(cell)[p].coord, 1e-13, "restart next position");
            check::vector(actual(cell)[p].velocity, expected(cell)[p].velocity, 1e-13, "restart next velocity");
        }
    }
    original->finalize();
    resumed->finalize();
}
}   // namespace

// Magnetized orbits with omega_c dt = 10 in a rough field: the E x B
// displacement responds to the field gradient with a factor ~ dt |grad E| / B
// near one, where plain Picard iteration cycles across cell planes and the
// solver falls back to substeps. The accelerated orbit iteration must reach
// the full-step Crank-Nicolson orbit for every particle.
void magnetized_orbits() {
    Domain domain;
    domain.init({8, 8, 8}, {1, 1, 1});
    Field3d electric(domain.size(), 3), magnetic(domain.size(), 3), current(domain.size(), 3);
    // A rough field indexed by the physical node, independent of the ghost layers.
    electric.setZero();
    magnetic.setZero();
    const int g = GHOST_CELLS;
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            for (int k = 0; k < 8; ++k)
                for (int c = 0; c < 3; ++c) {
                    electric(i + g, j + g, k + g, c) = 0.1 * std::sin(0.13 * (3 * ((i * 8 + j) * 8 + k) + c));
                    magnetic(i + g, j + g, k + g, c) = c == 2 ? 1.0 : 0.0;
                }
    algorithmsECCAPFIM::copy_periodic_nodes(electric, domain);
    algorithmsECCAPFIM::copy_periodic_nodes(magnetic, domain);
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> position(0.0, 8.0);
    std::normal_distribution<double> thermal(0.0, 0.1);
    const double dt = 10.0, qm = -1.0;
    double worst = 0.0;
    for (int n = 0; n < 500; ++n) {
        const Particle initial({position(rng), position(rng), position(rng)},
                               {thermal(rng), thermal(rng), thermal(rng)});
        current.setZero();
        const auto end = algorithmsECCAPFIM::advance_particle(initial, electric, magnetic, current, qm, -1e-3, dt,
                                                              domain, 1e-12, 100);
        Vector3R displacement = end.coord - initial.coord;
        for (int a = 0; a < 3; ++a) displacement[a] -= 8.0 * std::round(displacement[a] / 8.0);
        const auto fields =
            algorithmsECCAPFIM::gather_orbit(electric, magnetic, initial.coord, initial.coord + displacement, domain);
        const Vector3R vh = 0.5 * (initial.velocity + end.velocity);
        const double velocity_residual =
            (end.velocity - initial.velocity - dt * qm * (fields.electric + vh.cross(fields.magnetic))).length();
        worst = std::max({worst, velocity_residual, (displacement - dt * vh).length()});
    }
    check::require(worst < 1e-9,
                   "magnetized orbits reach the full-step CN orbit, worst residual " + std::to_string(worst));
}

int main() {
    return check::run([] {
        orbit_invariants();
        coupled_invariants();
        wall_orbits();
        accelerated_wall_orbits();
        decelerating_wall_orbits();
        mixed_boundary_deposit();
        wall_coupled_invariants();
        cold_plasma();
        check::scenario("ECCAPFIM magnetized orbits converge without substeps", magnetized_orbits);
        check::scenario("ECCAPFIM mass matrix equals frozen-orbit response", mass_matrix_response);
        check::scenario("ECCAPFIM mass-matrix preconditioners", [] { mass_matrix_preconditioners(Vector3I(6, 5, 7)); });
        check::scenario("ECCAPFIM mass-matrix preconditioners on a 1 x 1 x 8 column",
                        [] { mass_matrix_preconditioners(Vector3I(1, 1, 8)); });
        check::scenario("ECCAPFIM mass-matrix preconditioners with conducting z walls",
                        [] { mass_matrix_preconditioners(Vector3I(6, 5, 7), true); });
        maxwell_preconditioners();
        external_work();
        restart();
        particle_failure_dump();
        failure_is_not_a_step();
    });
}
