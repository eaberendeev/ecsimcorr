// Тесты общего механизма обрезки траектории частицы по граням домена:
// PathPolicy (какие грани удаляют частицы), разрешение политики из JSON
// граничных условий и shape_particle_path (остановка на первой пересечённой
// удаляющей грани). Сетка маленькая, без частиц.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "World.h"
#include "boundary_conditions.h"
#include "particle_path.h"
#include "vector3.h"

// Простой тестовый фреймворк без сторонних библиотек (как в domain.cpp).
namespace test {

static int passed = 0;
static int failed = 0;

void assert_true(bool condition, const std::string& test_name) {
    if (condition) {
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << std::endl;
        failed++;
    }
}

void assert_eq(double actual, double expected, double eps, const std::string& test_name) {
    if (std::abs(actual - expected) < eps) {
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << " (expected: " << expected << ", got: " << actual << ")" << std::endl;
        failed++;
    }
}

void assert_eq(int actual, int expected, const std::string& test_name) {
    if (actual == expected) {
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << " (expected: " << expected << ", got: " << actual << ")" << std::endl;
        failed++;
    }
}

void assert_face(Face actual, Face expected, const std::string& test_name) {
    assert_eq(static_cast<int>(actual), static_cast<int>(expected), test_name);
}

void print_summary() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

int get_failed() {
    return failed;
}

}   // namespace test

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
    test::assert_true(!policy.active(), "none(): inactive");
    test::assert_true(!policy.clips(Face::ZMIN), "none(): ZMIN not clipped");
    test::assert_true(!policy.clips(Face::UNAVAILABLE), "none(): UNAVAILABLE not clipped");

    policy.enable(Face::ZMIN);
    test::assert_true(policy.active(), "enable(ZMIN): active");
    test::assert_true(policy.clips(Face::ZMIN), "enable(ZMIN): ZMIN clipped");
    test::assert_true(!policy.clips(Face::ZMAX), "enable(ZMIN): ZMAX not clipped");
    test::assert_true(!policy.clips(Face::XMIN), "enable(ZMIN): XMIN not clipped");

    // Индекс UNAVAILABLE не должен приводить к выходу за пределы массива.
    test::assert_true(!policy.clips(Face::UNAVAILABLE), "enable(ZMIN): UNAVAILABLE not clipped");
    test::assert_true(!policy.clips(Face::Count), "enable(ZMIN): Count not clipped");
}

// =============================================================================
// Разрешение PathPolicy из JSON граничных условий
// =============================================================================
void test_handler_open() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json(R"([{"open": [{"face": "ZMIN"}, {"face": "ZMAX"}]}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    test::assert_true(policy.active(), "open ZMIN+ZMAX: active");
    test::assert_true(policy.clips(Face::ZMIN), "open ZMIN+ZMAX: clips ZMIN");
    test::assert_true(policy.clips(Face::ZMAX), "open ZMIN+ZMAX: clips ZMAX");
    test::assert_true(!policy.clips(Face::XMIN), "open ZMIN+ZMAX: does not clip XMIN");
    test::assert_true(!policy.clips(Face::XMAX), "open ZMIN+ZMAX: does not clip XMAX");
    test::assert_true(!policy.clips(Face::YMIN), "open ZMIN+ZMAX: does not clip YMIN");
    test::assert_true(!policy.clips(Face::YMAX), "open ZMIN+ZMAX: does not clip YMAX");

    BoundaryConditionHandler handler_x;
    handler_x.load_from_json(bc_json(R"([{"open": [{"face": "XMIN"}, {"face": "XMAX"}]}])"), domain);
    test::assert_true(handler_x.path_policy().clips(Face::XMIN), "open XMIN+XMAX: clips XMIN");
    test::assert_true(handler_x.path_policy().clips(Face::XMAX), "open XMIN+XMAX: clips XMAX");
}

void test_handler_periodic() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json(R"([{"periodic": {"face": "XMIN"}}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    test::assert_true(!policy.clips(Face::XMIN), "periodic X: does not clip XMIN");
    test::assert_true(!policy.clips(Face::XMAX), "periodic X: does not clip XMAX");
    test::assert_true(!policy.active(), "periodic only: inactive");
}

void test_handler_reload() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    const auto periodic = bc_json(R"([{"periodic": {"face": "XMIN"}}])");
    handler.load_from_json(periodic, domain);
    handler.load_from_json(bc_json("[]"), domain);
    test::assert_true(!handler.is_periodic(0), "reload clears periodic axis");
    handler.load_from_json(periodic, domain);
    Field3d first(domain.size(), 1), second(domain.size(), 1);
    first.data().setOnes();
    second.data().setOnes();
    handler.apply_to_fields(first, FieldType::DENSITY, domain);
    handler.load_from_json(periodic, domain);
    handler.apply_to_fields(second, FieldType::DENSITY, domain);
    test::assert_eq((first - second).norm(), 0.0, kEps, "reload retains periodic field folding");
    handler.load_from_json(nlohmann::json::object(), domain);
    test::assert_true(!handler.is_periodic(0), "absent boundary list clears periodic axis");
}

void test_handler_absent() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json("[]"), domain);
    test::assert_true(!handler.path_policy().active(), "no boundary conditions: inactive");
    test::assert_true(!handler.path_policy().clips(Face::ZMIN), "no boundary conditions: ZMIN not clipped");
}

void test_handler_er0() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json(R"([{"er0": {"face": "ZMIN"}}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    test::assert_true(policy.active(), "er0 ZMIN: active");
    test::assert_true(policy.clips(Face::ZMIN), "er0 ZMIN: clips ZMIN (inherits open removal)");
    test::assert_true(!policy.clips(Face::ZMAX), "er0 ZMIN: does not clip ZMAX");
}

void test_handler_bphi() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(
        bc_json(R"([{"bphi": {"face": "ZMAX", "gap": 0.0, "radius": 1.0, "electron_threshold_energy_": 0.0}}])"),
        domain);
    const PathPolicy& policy = handler.path_policy();
    test::assert_true(policy.active(), "bphi ZMAX: active");
    test::assert_true(policy.clips(Face::ZMAX), "bphi ZMAX: clips ZMAX");
    test::assert_true(!policy.clips(Face::ZMIN), "bphi ZMAX: does not clip ZMIN");
}

void test_handler_electron_reflection() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json(R"([{"electron_reflection": {"face": "ZMIN"}}])"), domain);
    const PathPolicy& policy = handler.path_policy();
    test::assert_true(policy.active(), "electron_reflection ZMIN: active");
    test::assert_true(policy.clips(Face::ZMIN), "electron_reflection ZMIN: clips ZMIN");
    test::assert_true(!policy.clips(Face::ZMAX), "electron_reflection ZMIN: does not clip ZMAX");
}

// =============================================================================
// shape_particle_path: неактивная политика
// =============================================================================
void test_shape_inactive() {
    Domain domain = make_domain();
    const Vector3R start(2.0, 2.0, 1.5);
    const Vector3R trial_end(2.0, 2.0, -0.5);

    Vector3R end_for_current(0.0, 0.0, 0.0);
    const PathResult result = shape_particle_path(start, trial_end, domain, PathPolicy::none(), end_for_current);
    test::assert_true(!result.clipped, "inactive: not clipped");
    test::assert_face(result.face, Face::UNAVAILABLE, "inactive: face UNAVAILABLE");
    test::assert_eq(end_for_current.x(), trial_end.x(), kEps, "inactive: endpoint x unchanged");
    test::assert_eq(end_for_current.y(), trial_end.y(), kEps, "inactive: endpoint y unchanged");
    test::assert_eq(end_for_current.z(), trial_end.z(), kEps, "inactive: endpoint z unchanged");
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
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(result.clipped, "ZMIN: clipped");
    test::assert_face(result.face, Face::ZMIN, "ZMIN: face");
    test::assert_eq(end_for_current.z(), domain.geom.box_min.z(), kEps, "ZMIN: z snapped to box_min");
    test::assert_eq(end_for_current.x(), start.x(), kEps, "ZMIN: x unchanged");
    test::assert_eq(end_for_current.y(), start.y(), kEps, "ZMIN: y unchanged");
}

void test_shape_zmax() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMAX);

    const Vector3R start(2.0, 2.0, 3.5);
    const Vector3R trial_end(2.0, 2.0, 4.5);
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(result.clipped, "ZMAX: clipped");
    test::assert_face(result.face, Face::ZMAX, "ZMAX: face");
    test::assert_eq(end_for_current.z(), domain.geom.box_max.z(), kEps, "ZMAX: z snapped to box_max");
}

void test_shape_xmin_only() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::XMIN);

    const Vector3R start(1.5, 2.0, 2.0);
    const Vector3R trial_end(-0.5, 2.0, 2.0);
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(result.clipped, "XMIN: clipped");
    test::assert_face(result.face, Face::XMIN, "XMIN: face");
    test::assert_eq(end_for_current.x(), domain.geom.box_min.x(), kEps, "XMIN: x snapped to box_min");
}

void test_shape_xmax_not_enabled() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::XMAX);

    // Тот же отрезок, что и в XMIN-тесте: он пересекает XMIN, а не XMAX.
    const Vector3R start(1.5, 2.0, 2.0);
    const Vector3R trial_end(-0.5, 2.0, 2.0);
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(!result.clipped, "XMAX-only: XMIN crossing not clipped");
    test::assert_face(result.face, Face::UNAVAILABLE, "XMAX-only: face UNAVAILABLE");
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
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(result.clipped, "corner: clipped");
    test::assert_face(result.face, Face::XMAX, "corner: first crossed face is XMAX");
    test::assert_eq(end_for_current.x(), domain.geom.box_max.x(), kEps, "corner: x snapped to box_max");
    test::assert_true(end_for_current.z() < domain.geom.box_max.z(), "corner: z stopped before ZMAX");
}

void test_shape_start_on_plane_inward() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMIN);

    const Vector3R start(2.0, 2.0, domain.geom.box_min.z());
    const Vector3R trial_end(2.0, 2.0, 1.0);
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(!result.clipped, "start on ZMIN moving inward: not clipped");
    test::assert_eq(end_for_current.z(), trial_end.z(), kEps, "start on ZMIN moving inward: endpoint unchanged");
}

void test_shape_parallel_to_wall() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMIN);

    const Vector3R start(1.0, 1.0, domain.geom.box_min.z());
    const Vector3R trial_end(2.0, 2.0, domain.geom.box_min.z());
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(!result.clipped, "segment parallel to ZMIN: not clipped");
    test::assert_eq(end_for_current.z(), trial_end.z(), kEps, "segment parallel to ZMIN: endpoint unchanged");
}

void test_shape_incoming_from_outside() {
    Domain domain = make_domain();
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::ZMIN);

    const Vector3R start(2.0, 2.0, -1.0);
    const Vector3R trial_end(2.0, 2.0, 1.0);
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(!result.clipped, "incoming from outside ZMIN: not clipped");
    test::assert_eq(end_for_current.z(), trial_end.z(), kEps, "incoming from outside: endpoint unchanged");
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
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(result.clipped, "cylinder: clipped");
    test::assert_face(result.face, Face::CYLINDER, "cylinder: face");
    const double distance = std::hypot(end_for_current.x() - domain.geom.cyl_center.x(),
                                       end_for_current.y() - domain.geom.cyl_center.y());
    test::assert_eq(distance, domain.geom.cyl_radius, kEps, "cylinder: endpoint on |p-c| = R");
    test::assert_eq(end_for_current.z(), start.z(), kEps, "cylinder: z unchanged");
}

void test_shape_cylinder_inside() {
    Domain domain = make_domain();
    domain.set_cylinder(Vector3R(0.0, 0.0, 0.0), 1.0);
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::CYLINDER);

    const Vector3R start(0.0, 0.0, 2.0);
    const Vector3R trial_end(0.5, 0.0, 2.0);
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(!result.clipped, "cylinder: segment staying inside not clipped");
    test::assert_eq(end_for_current.x(), trial_end.x(), kEps, "cylinder inside: endpoint x unchanged");
}

void test_shape_cylinder_without_geometry() {
    // Политика включает CYLINDER, но в геометрии цилиндра нет.
    Domain domain = make_domain();
    test::assert_true(!domain.geom.use_cylinder, "no-cylinder domain: use_cylinder false");
    PathPolicy policy = PathPolicy::none();
    policy.enable(Face::CYLINDER);

    const Vector3R start(0.0, 0.0, 2.0);
    const Vector3R trial_end(2.0, 0.0, 2.0);
    Vector3R end_for_current = trial_end;
    const PathResult result = shape_particle_path(start, trial_end, domain, policy, end_for_current);

    test::assert_true(!result.clipped, "CYLINDER policy without cylinder geometry: not clipped");
    test::assert_face(result.face, Face::UNAVAILABLE, "CYLINDER policy without geometry: face UNAVAILABLE");
    test::assert_eq(end_for_current.x(), trial_end.x(), kEps, "CYLINDER policy without geometry: endpoint unchanged");
}

}   // namespace

int main() {
    // PathPolicy
    test_path_policy_accessors();

    // Разрешение политики из JSON
    test_handler_open();
    test_handler_periodic();
    test_handler_reload();
    test_handler_absent();
    test_handler_er0();
    test_handler_bphi();
    test_handler_electron_reflection();

    // shape_particle_path: плоские грани
    test_shape_inactive();
    test_shape_zmin();
    test_shape_zmax();
    test_shape_xmin_only();
    test_shape_xmax_not_enabled();
    test_shape_corner_first_face();
    test_shape_start_on_plane_inward();
    test_shape_parallel_to_wall();
    test_shape_incoming_from_outside();

    // shape_particle_path: цилиндр
    test_shape_cylinder_exit();
    test_shape_cylinder_inside();
    test_shape_cylinder_without_geometry();

    test::print_summary();

    return test::get_failed() == 0 ? 0 : 1;
}
