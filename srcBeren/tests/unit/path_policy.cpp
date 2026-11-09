#include "assert_contract.h"
// Тесты общего механизма обрезки траектории частицы по граням домена:
// PathPolicy (какие грани удаляют частицы), разрешение политики из JSON
// граничных условий и shape_particle_path (остановка на первой пересечённой
// удаляющей грани), а также последующего размещения частиц по ячейкам.
#include <cmath>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "config/boundary_config.h"
#include "core/vector3.h"
#include "grid/domain.h"
#include "particles/particle_path.h"
#include "particles/particle_sources.h"
#include "simulation/boundary_conditions.h"
#include "test_support.h"

namespace {

constexpr double kEps = 1e-12;

// Домен 4x4x4 ячеек по 1.0: box [0,4]^3.
Domain make_domain() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1.0, 1.0, 1.0));
    return domain;
}

nlohmann::json bc_json(const std::string& conditions) {
    return nlohmann::json::parse(std::string(R"({"Boundary_conditions": )") + conditions + "}");
}

// =============================================================================
// PathPolicy
// =============================================================================
void test_path_policy_accessors() {
    PathPolicy policy = PathPolicy::none();
    check::require(!policy.active(), "none(): inactive");
    check::require(!policy.clips(Face::ZMIN), "none(): ZMIN not clipped");
    check::require(!policy.clips(Face::UNAVAILABLE), "none(): UNAVAILABLE not clipped");

    policy.enable(Face::ZMIN);
    check::require(policy.active(), "enable(ZMIN): active");
    check::require(policy.clips(Face::ZMIN), "enable(ZMIN): ZMIN clipped");
    check::require(!policy.clips(Face::ZMAX), "enable(ZMIN): ZMAX not clipped");
    check::require(!policy.clips(Face::XMIN), "enable(ZMIN): XMIN not clipped");

    // Индекс UNAVAILABLE не должен приводить к выходу за пределы массива.
    check::require(!policy.clips(Face::UNAVAILABLE), "enable(ZMIN): UNAVAILABLE not clipped");
    check::require(!policy.clips(Face::Count), "enable(ZMIN): Count not clipped");
}

// =============================================================================
// Разрешение PathPolicy из JSON граничных условий
// =============================================================================
void test_handler_open() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(bc_json(R"([{"open": [{"face": "ZMIN"}, {"face": "ZMAX"}]}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    check::require(policy.active(), "open ZMIN+ZMAX: active");
    check::require(policy.clips(Face::ZMIN), "open ZMIN+ZMAX: clips ZMIN");
    check::require(policy.clips(Face::ZMAX), "open ZMIN+ZMAX: clips ZMAX");
    check::require(!policy.clips(Face::XMIN), "open ZMIN+ZMAX: does not clip XMIN");
    check::require(!policy.clips(Face::XMAX), "open ZMIN+ZMAX: does not clip XMAX");
    check::require(!policy.clips(Face::YMIN), "open ZMIN+ZMAX: does not clip YMIN");
    check::require(!policy.clips(Face::YMAX), "open ZMIN+ZMAX: does not clip YMAX");

    BoundaryConditionHandler handler_x;
    handler_x = boundary_conditions_from_config(bc_json(R"([{"open": [{"face": "XMIN"}, {"face": "XMAX"}]}])"), domain);
    check::require(handler_x.path_policy().clips(Face::XMIN), "open XMIN+XMAX: clips XMIN");
    check::require(handler_x.path_policy().clips(Face::XMAX), "open XMIN+XMAX: clips XMAX");
}

void test_handler_periodic() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(bc_json(R"([{"periodic": {"face": "XMIN"}}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    check::require(!policy.clips(Face::XMIN), "periodic X: does not clip XMIN");
    check::require(!policy.clips(Face::XMAX), "periodic X: does not clip XMAX");
    check::require(!policy.active(), "periodic only: inactive");
}

void test_handler_reload() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    const auto periodic = bc_json(R"([{"periodic": {"face": "XMIN"}}])");
    handler = boundary_conditions_from_config(periodic, domain);
    handler = boundary_conditions_from_config(bc_json("[]"), domain);
    check::require(!handler.is_periodic(0), "reload clears periodic axis");
    handler = boundary_conditions_from_config(periodic, domain);
    Field3d first(domain.size(), 1), second(domain.size(), 1);
    first.data().setOnes();
    second.data().setOnes();
    handler.apply_to_fields(first, FieldType::DENSITY, domain);
    handler = boundary_conditions_from_config(periodic, domain);
    handler.apply_to_fields(second, FieldType::DENSITY, domain);
    check::equal((first - second).norm(), 0.0, kEps, "reload retains periodic field folding");
    handler = boundary_conditions_from_config(nlohmann::json::object(), domain);
    check::require(!handler.is_periodic(0), "absent boundary list clears periodic axis");
}

void test_handler_absent() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(bc_json("[]"), domain);
    check::require(!handler.path_policy().active(), "no boundary conditions: inactive");
    check::require(!handler.path_policy().clips(Face::ZMIN), "no boundary conditions: ZMIN not clipped");
}

void test_handler_er0() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(bc_json(R"([{"er0": {"face": "ZMIN"}}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    check::require(policy.active(), "er0 ZMIN: active");
    check::require(policy.clips(Face::ZMIN), "er0 ZMIN: clips ZMIN (inherits open removal)");
    check::require(!policy.clips(Face::ZMAX), "er0 ZMIN: does not clip ZMAX");
}

void test_handler_bphi() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(
        bc_json(R"([{"bphi": {"face": "ZMAX", "gap": 0.0, "radius": 1.0, "electron_threshold_energy_": 0.0}}])"),
        domain);
    const PathPolicy& policy = handler.path_policy();
    check::require(policy.active(), "bphi ZMAX: active");
    check::require(policy.clips(Face::ZMAX), "bphi ZMAX: clips ZMAX");
    check::require(!policy.clips(Face::ZMIN), "bphi ZMAX: does not clip ZMIN");
}

void test_handler_electron_reflection() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(bc_json(R"([{"electron_reflection": {"face": "ZMIN"}}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    check::require(policy.active(), "electron_reflection ZMIN: active");
    check::require(policy.clips(Face::ZMIN), "electron_reflection ZMIN: clips ZMIN");
    check::require(!policy.clips(Face::ZMAX), "electron_reflection ZMIN: does not clip ZMAX");
}

void test_selective_wall_law() {
    const Domain domain = make_domain();
    for (Face face : {Face::ZMIN, Face::ZMAX}) {
        const std::string face_name = face == Face::ZMIN ? "ZMIN" : "ZMAX";
        BoundaryConditionHandler bphi = boundary_conditions_from_config(
            {{"Boundary_conditions",
              nlohmann::json::array(
                  {{{"bphi",
                     {{"face", face_name}, {"gap", 0}, {"radius", 1}, {"electron_threshold_energy_", 0.125}}}}})}},
            domain);
        ElectronReflectionCondition reflection(face, 1, 0.125);
        for (const std::string name : {"Electrons", "Ions"}) {
            ParticlesArray species(check::species(name, name == "Electrons" ? -1 : 1, 1), domain);
            for (double x : {2.0, 3.0, 3.01}) {
                for (double speed : {0.2, 0.5, 0.6}) {
                    const double z = face == Face::ZMIN ? -0.01 : 4.01;
                    const Particle particle(Vector3R(x, 2, z), Vector3R(0, 0, speed));
                    // Radius and energy thresholds include their endpoints:
                    // m=1, Kz=vz^2/2 <= 0.125, radius <= 1 about (2,2).
                    const ParticleFate expected =
                        name == "Electrons" && x <= 3 && speed <= 0.5 ? ParticleFate::Reflected : ParticleFate::Removed;
                    for (Face selected : {Face::UNAVAILABLE, face}) {
                        const ParticleFateResult a = bphi.classify_particle(particle, species, domain, selected);
                        const ParticleFateResult b = reflection.classify_particle(particle, species, domain, selected);
                        check::require(a.fate == expected && a.face == face, "bphi radius/energy/species wall law");
                        check::require(b.fate == expected && b.face == face, "electron reflection wall law");
                    }
                }
            }
        }
    }
}

// =============================================================================
// shape_particle_path: неактивная политика
// =============================================================================
void test_shape_inactive() {
    Domain domain = make_domain();
    const Vector3R start(2.0, 2.0, 1.5);
    const Vector3R trial_end(2.0, 2.0, -0.5);

    Vector3R end_for_current(0.0, 0.0, 0.0);
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, PathPolicy::none(), end_for_current);
    check::require(!result.clipped, "inactive: not clipped");
    check::equal(result.face, Face::UNAVAILABLE, "inactive: face UNAVAILABLE");
    check::equal(end_for_current.x(), trial_end.x(), kEps, "inactive: endpoint x unchanged");
    check::equal(end_for_current.y(), trial_end.y(), kEps, "inactive: endpoint y unchanged");
    check::equal(end_for_current.z(), trial_end.z(), kEps, "inactive: endpoint z unchanged");
}

// =============================================================================
// shape_particle_path: плоские грани
// =============================================================================
void test_shape_zmin() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMIN);

    const Vector3R start(2.0, 2.0, 1.5);
    const Vector3R trial_end(2.0, 2.0, -0.5);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(result.clipped, "ZMIN: clipped");
    check::equal(result.face, Face::ZMIN, "ZMIN: face");
    check::equal(end_for_current.z(), domain.geom.box_min.z(), kEps, "ZMIN: z snapped to box_min");
    check::equal(end_for_current.x(), start.x(), kEps, "ZMIN: x unchanged");
    check::equal(end_for_current.y(), start.y(), kEps, "ZMIN: y unchanged");
}

void test_shape_zmax() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMAX);

    const Vector3R start(2.0, 2.0, 3.5);
    const Vector3R trial_end(2.0, 2.0, 4.5);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(result.clipped, "ZMAX: clipped");
    check::equal(result.face, Face::ZMAX, "ZMAX: face");
    check::equal(end_for_current.z(), domain.geom.box_max.z(), kEps, "ZMAX: z snapped to box_max");
}

void test_shape_xmin_only() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::XMIN);

    const Vector3R start(1.5, 2.0, 2.0);
    const Vector3R trial_end(-0.5, 2.0, 2.0);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(result.clipped, "XMIN: clipped");
    check::equal(result.face, Face::XMIN, "XMIN: face");
    check::equal(end_for_current.x(), domain.geom.box_min.x(), kEps, "XMIN: x snapped to box_min");
}

void test_shape_xmax_not_enabled() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::XMAX);

    // Тот же отрезок, что и в XMIN-тесте: он пересекает XMIN, а не XMAX.
    const Vector3R start(1.5, 2.0, 2.0);
    const Vector3R trial_end(-0.5, 2.0, 2.0);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(!result.clipped, "XMAX-only: XMIN crossing not clipped");
    check::equal(result.face, Face::UNAVAILABLE, "XMAX-only: face UNAVAILABLE");
}

void test_shape_corner_first_face() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::XMAX);
    policy.enable(Face::ZMAX);

    // XMAX достигается при t = 1/3 (x: 3.5 -> 4), ZMAX позже (t = 0.5).
    const Vector3R start(3.5, 2.0, 3.5);
    const Vector3R trial_end(5.0, 2.0, 4.5);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(result.clipped, "corner: clipped");
    check::equal(result.face, Face::XMAX, "corner: first crossed face is XMAX");
    check::equal(end_for_current.x(), domain.geom.box_max.x(), kEps, "corner: x snapped to box_max");
    check::require(end_for_current.z() < domain.geom.box_max.z(), "corner: z stopped before ZMAX");
}

void test_start_on_each_wall() {
    const auto domain = make_domain();
    const Face lower[] = {Face::XMIN, Face::YMIN, Face::ZMIN};
    const Face upper[] = {Face::XMAX, Face::YMAX, Face::ZMAX};
    for (int axis = 0; axis < 3; ++axis) {
        for (bool high : {false, true}) {
            const Face face = high ? upper[axis] : lower[axis];
            PathPolicy policy;
            policy.enable(face);
            Vector3R start(2, 2, 2), end;
            start[axis] = high ? 4 : 0;
            Vector3R outward = start;
            outward[axis] += high ? 0.2 : -0.2;
            const auto hit = shape_particle_path(start, outward, domain, policy, end);
            check::require(hit.clipped && hit.face == face, "outward start must hit its wall");
            check::near(hit.t_hit, 0, 0, "wall start has t_hit=0");
            check::vector(end, start, 0, "wall start has zero internal trajectory");
            const auto inward = start - (outward - start);
            check::require(!shape_particle_path(start, inward, domain, policy, end).clipped,
                           "inward movement from wall must remain free");
        }
    }
}

void test_shape_parallel_to_wall() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMIN);

    const Vector3R start(1.0, 1.0, domain.geom.box_min.z());
    const Vector3R trial_end(2.0, 2.0, domain.geom.box_min.z());
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(!result.clipped, "segment parallel to ZMIN: not clipped");
    check::equal(end_for_current.z(), trial_end.z(), kEps, "segment parallel to ZMIN: endpoint unchanged");
}

void test_shape_incoming_from_outside() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMIN);

    const Vector3R start(2.0, 2.0, -1.0);
    const Vector3R trial_end(2.0, 2.0, 1.0);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(!result.clipped, "incoming from outside ZMIN: not clipped");
    check::equal(end_for_current.z(), trial_end.z(), kEps, "incoming from outside: endpoint unchanged");
}

// =============================================================================
// shape_particle_path: цилиндр
// =============================================================================
void test_shape_cylinder_exit() {
    Domain domain = make_domain();
    domain.set_cylinder(Vector3R(0.0, 0.0, 0.0), 1.0);
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::CYLINDER);

    const Vector3R start(0.0, 0.0, 2.0);
    const Vector3R trial_end(2.0, 0.0, 2.0);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(result.clipped, "cylinder: clipped");
    check::equal(result.face, Face::CYLINDER, "cylinder: face");
    const double distance =
        std::hypot(end_for_current.x() - domain.geom.cyl_center.x(), end_for_current.y() - domain.geom.cyl_center.y());
    check::equal(distance, domain.geom.cyl_radius, kEps, "cylinder: endpoint on |p-c| = R");
    check::equal(end_for_current.z(), start.z(), kEps, "cylinder: z unchanged");
}

void test_shape_cylinder_inside() {
    Domain domain = make_domain();
    domain.set_cylinder(Vector3R(0.0, 0.0, 0.0), 1.0);
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::CYLINDER);

    const Vector3R start(0.0, 0.0, 2.0);
    const Vector3R trial_end(0.5, 0.0, 2.0);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(!result.clipped, "cylinder: segment staying inside not clipped");
    check::equal(end_for_current.x(), trial_end.x(), kEps, "cylinder inside: endpoint x unchanged");
}

void test_shape_cylinder_without_geometry() {
    // Политика включает CYLINDER, но в геометрии цилиндра нет.
    Domain domain = make_domain();
    check::require(!domain.geom.use_cylinder, "no-cylinder domain: use_cylinder false");
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::CYLINDER);

    const Vector3R start(0.0, 0.0, 2.0);
    const Vector3R trial_end(2.0, 0.0, 2.0);
    Vector3R end_for_current = trial_end;
    const BoundaryCrossing result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    check::require(!result.clipped, "CYLINDER policy without cylinder geometry: not clipped");
    check::equal(result.face, Face::UNAVAILABLE, "CYLINDER policy without geometry: face UNAVAILABLE");
    check::equal(end_for_current.x(), trial_end.x(), kEps, "CYLINDER policy without geometry: endpoint unchanged");
}

void particle_binning_checks_storage_before_access() {
    {
        constexpr int ghost = GHOST_CELLS;
        Domain domain;
        domain.init(Vector3I(3, 4, 5), Vector3R(0.5, 1, 2));
        ParticlesArray particles(check::species("Electrons"), domain);
        for (int axis = 0; axis < 3; ++axis) {
            for (double cell : {-ghost - 0.125, double(domain.num_cells(axis) + ghost + 1)}) {
                Vector3R position(0.25, 0.5, 1);
                position[axis] = cell * domain.cell_size(axis);
                check::assertion_fails([&] { particles.add_particle(Particle(position, Vector3R(0.0))); });
            }
        }
        check::equal(particles.get_total_num_of_particles(), 0, "invalid insertions must not change particle count");
    }
    {
        constexpr int ghost = GHOST_CELLS;
        Domain domain;
        domain.init(Vector3I(3, 4, 5), Vector3R(0.5, 1, 2));
        ParticlesArray particles(check::species("Electrons"), domain);
        const Vector3I start(1, ghost + 1, ghost + 1);
        particles.add_particle(Particle(Vector3R((1.25 - ghost) * 0.5, 1.5, 3), Vector3R(0.0)));
        // A neighboring-bin move reaches raw zero for every ghost width,
        // including the lower halo when present.
        particles(start.x(), start.y(), start.z()).front().coord.x() = (0.25 - ghost) * 0.5;
        particles.update_cells(domain);
        check::equal(particles(0, ghost + 1, ghost + 1).size(), size_t(1),
                     "rebinning into lower allocated cell preserves the particle");
        check::equal(particles.get_total_num_of_particles(), 1, "neighbor rebinning conserves particle count");
    }
}

struct FixedSourceVelocity : IVelocityDistribution {
    Vector3R sample(LehmerEngine&) const override {
        return Vector3R(0.0);
    }

    Vector3R sample(const Vector3R&, LehmerEngine&) const override {
        return Vector3R(0.0);
    }
};

void particle_sources_check_serial_and_parallel_storage() {
    Domain domain;
    domain.init(Vector3I(3, 4, 5), Vector3R(0.5, 1, 2));
    for (int count : {1, 15'001}) {
        check::assertion_fails([&] {
            ParticlesArray particles(check::species("Electrons"), domain);
            InitialDistribution source(
                std::make_shared<RectangleDistribution>(Vector3R((-GHOST_CELLS - 0.25) * 0.5, 1, 2), Vector3R(0.0)),
                std::make_shared<FixedSourceVelocity>(), "initial", count, 1, 1);
            LehmerEngine space(13), velocity(15);
            add_particles_from_distribution(particles, source, space, velocity, domain, 0, false);
        });
    }
}

}   // namespace

int main() {
    return check::run([] {
        // Fork assertion checks before this process starts any OpenMP workers.
        check::scenario("Serial and parallel particle sources reject positions outside allocated storage",
                        particle_sources_check_serial_and_parallel_storage);
        check::scenario(
            "Particle insertion and neighboring-bin motion respect halo bounds for the configured ghost width",
            particle_binning_checks_storage_before_access);

        check::scenario("Enabled absorbing faces and invalid face indices", test_path_policy_accessors);
        check::scenario("Open boundary settings clip only the configured faces", test_handler_open);
        check::scenario("Periodic faces do not terminate the current path", test_handler_periodic);
        check::scenario("Reloading boundary settings clears old state and folds periodic fields once",
                        test_handler_reload);
        check::scenario("No boundary conditions means no path clipping", test_handler_absent);
        check::scenario("Er0 retains its absorbing particle boundary", test_handler_er0);
        check::scenario("Bphi contributes the configured particle boundary", test_handler_bphi);
        check::scenario("Selective electron wall appears in the path policy", test_handler_electron_reflection);
        check::scenario("Both selective walls obey the radius, energy and species criteria", test_selective_wall_law);
        check::scenario("Inactive path policy preserves the trial trajectory", test_shape_inactive);
        check::scenario("Absorbed path ends at ZMIN", test_shape_zmin);
        check::scenario("Absorbed path ends at ZMAX", test_shape_zmax);
        check::scenario("Absorbed path ends at XMIN", test_shape_xmin_only);
        check::scenario("An unconfigured wall cannot clip a path", test_shape_xmax_not_enabled);
        check::scenario("The first geometrically crossed face determines the endpoint", test_shape_corner_first_face);
        check::scenario("Starting on a wall: outward hits at t=0, inward remains free", test_start_on_each_wall);
        check::scenario("Tangential movement along a wall is not an outgoing crossing", test_shape_parallel_to_wall);
        check::scenario("An incoming path is not clipped as an outgoing path", test_shape_incoming_from_outside);
        check::scenario("Cylinder crossing ends at the known radial surface", test_shape_cylinder_exit);
        check::scenario("A path inside the cylinder is unchanged", test_shape_cylinder_inside);
        check::scenario("A cylindrical flag without cylindrical geometry cannot clip a path",
                        test_shape_cylinder_without_geometry);
    });
}
