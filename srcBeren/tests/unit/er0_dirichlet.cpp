// Тесты электрической Дирихле-маски условия er0 (режим "solv") и утилит
// наложения строк матрицы / правой части. Сетка маленькая, без частиц.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "sgs.h"
#include "types.h"

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

constexpr double kEps = 1e-15;
// Тот же допуск, что у OpenBoundaryCondition::eps_ (см. boundary_conditions.h).
constexpr double kMaskEps = 1e-10;

// Домен 8x8x6 ячеек по 1.0, ghost=1 -> сетка 11x11x9.
// Плоскость Ex/Ey: z = (k - ghost) * dz; плоскость Ez: z = (k - ghost + 0.5) * dz.
Domain make_domain(bool with_cylinder = false) {
    Domain domain;
    domain.init(Vector3I(8, 8, 6), Vector3R(1.0, 1.0, 1.0));
    if (with_cylinder) {
        domain.set_cylinder(Vector3R(4.0, 4.0, 0.0), 4.0);
    }
    return domain;
}

nlohmann::json bc_json(const std::string& conditions) {
    return nlohmann::json::parse(std::string(R"({"Boundary_conditions": )") + conditions + "}");
}

struct Er0Params {
    double inner_radius = 1.0;
    double width = 1.0;
    double potential_drop = 0.1 / SGS::MC2;
};

int face_plane_k(const Domain& domain, Face face) {
    const int ghost = domain.grid.ghost_cells();
    return (face == Face::ZMIN) ? ghost : domain.size().z() - ghost - 1;
}

// Независимый расчёт предписанного значения в узле (i,j,k,d).
double expected_value(const Domain& domain, Face face, int i, int j, int k, int d, const Er0Params& p) {
    if (d >= 2 || k != face_plane_k(domain, face)) {
        return 0.0;
    }
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double center_x = 0.5 * domain.num_cells().x() * dx;
    const double center_y = 0.5 * domain.num_cells().y() * dy;
    const Vector3R pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
    const double xx = pos.x() - center_x;
    const double yy = pos.y() - center_y;
    const double rr = std::hypot(xx, yy);
    if (rr < p.inner_radius || rr > p.inner_radius + p.width || rr <= 0.0) {
        return 0.0;
    }
    return -p.potential_drop / p.width * ((d == 0) ? xx : yy) / rr;
}

bool in_mask(const Domain& domain, Face face, int i, int j, int k, int d) {
    const Vector3R pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
    return domain.geom.is_outside_face(face, pos, -0.0 + kMaskEps) &&
           domain.geom.contains_ignoring_face(face, pos, 0.0);
}

struct ExpectedMask {
    std::vector<ElectricDirichletEntry> entries;
    int ring_count = 0;
};

// Независимая реконструкция ожидаемой маски по геометрическим предикатам.
ExpectedMask build_expected(const Domain& domain, Face face, const Er0Params& p) {
    ExpectedMask result;
    for (int i = 0; i < domain.size().x(); ++i) {
        for (int j = 0; j < domain.size().y(); ++j) {
            for (int k = 0; k < domain.size().z(); ++k) {
                for (int d = 0; d < 3; ++d) {
                    if (!in_mask(domain, face, i, j, k, d)) {
                        continue;
                    }
                    const double value = expected_value(domain, face, i, j, k, d, p);
                    result.entries.push_back({domain.vind(i, j, k, d), value});
                    if (value != 0.0) {
                        result.ring_count++;
                    }
                }
            }
        }
    }
    std::sort(result.entries.begin(), result.entries.end(),
              [](const ElectricDirichletEntry& a, const ElectricDirichletEntry& b) { return a.row < b.row; });
    return result;
}

void check_mask_matches(const std::vector<ElectricDirichletEntry>& actual, const ExpectedMask& expected,
                        const std::string& prefix) {
    test::assert_eq(static_cast<int>(actual.size()), static_cast<int>(expected.entries.size()),
                    prefix + ": размер маски");
    if (actual.size() != expected.entries.size()) {
        return;
    }
    int row_mismatch = 0;
    int value_mismatch = 0;
    for (size_t n = 0; n < actual.size(); ++n) {
        if (actual[n].row != expected.entries[n].row) {
            row_mismatch++;
        }
        if (std::abs(actual[n].value - expected.entries[n].value) > kEps) {
            value_mismatch++;
        }
    }
    test::assert_eq(row_mismatch, 0, prefix + ": набор строк");
    test::assert_eq(value_mismatch, 0, prefix + ": значения");

    int nonzero = 0;
    for (const auto& e : actual) {
        if (e.value != 0.0) {
            nonzero++;
        }
    }
    test::assert_eq(nonzero, expected.ring_count, prefix + ": число ненулевых (кольцевых) DOF");
}

void test_mask_zmin() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json(R"([{"open": [{"face": "ZMAX"}]}, {"er0": {"face": "ZMIN", )"
                                   R"("inner_radius": 1.0, "width": 1.0, "potential_drop": 0.1, "mode": "solv"}}])"),
                           domain);
    const auto& entries = handler.electric_dirichlet(domain);
    const Er0Params p;
    check_mask_matches(entries, build_expected(domain, Face::ZMIN, p), "ZMIN");

    // «Лишний» узел k=0 (z<0) входит в маску для всех компонент и равен нулю.
    int k0_count = 0;
    int k0_nonzero = 0;
    for (const auto& e : entries) {
        const auto [i, j, k, d] = domain.grid.pos_vind_range(e.row);
        (void) i;
        (void) j;
        (void) d;
        if (k == 0) {
            k0_count++;
            if (e.value != 0.0) {
                k0_nonzero++;
            }
        }
    }
    test::assert_true(k0_count > 0, "ZMIN: лишняя плоскость k=0 в маске");
    test::assert_eq(k0_nonzero, 0, "ZMIN: лишняя плоскость k=0 обнулена");

    // Ez на k=1 (z=+0.5) лежит внутри домена — это свободный DOF, не в маске.
    int ez_inside = 0;
    for (int i = 0; i < domain.size().x(); ++i) {
        for (int j = 0; j < domain.size().y(); ++j) {
            if (in_mask(domain, Face::ZMIN, i, j, 1, 2)) {
                ez_inside++;
            }
        }
    }
    test::assert_eq(ez_inside, 0, "ZMIN: Ez(k=1) не в маске");
}

void test_mask_zmax() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json(R"([{"open": [{"face": "ZMIN"}]}, {"er0": {"face": "ZMAX", )"
                                   R"("inner_radius": 1.0, "width": 1.0, "potential_drop": 0.1, "mode": "solv"}}])"),
                           domain);
    const auto& entries = handler.electric_dirichlet(domain);
    const Er0Params p;
    check_mask_matches(entries, build_expected(domain, Face::ZMAX, p), "ZMAX");

    const int k_face = face_plane_k(domain, Face::ZMAX);
    const int k_extra = domain.size().z() - 1;

    int ring_face = 0;
    int extra = 0;
    int ez_face = 0;
    int zero_violations = 0;
    for (const auto& e : entries) {
        const auto [i, j, k, d] = domain.grid.pos_vind_range(e.row);
        (void) i;
        (void) j;
        if (k == k_face && d < 2 && e.value != 0.0) {
            ring_face++;
        }
        if (k == k_face && d == 2) {
            ez_face++;
            if (e.value != 0.0) {
                zero_violations++;
            }
        }
        if (k == k_extra) {
            extra++;
            if (e.value != 0.0) {
                zero_violations++;
            }
        }
    }
    test::assert_true(ring_face > 0, "ZMAX: кольцо на k=size_z-2");
    test::assert_true(extra > 0, "ZMAX: лишняя плоскость k=size_z-1 в маске");
    test::assert_true(ez_face > 0, "ZMAX: Ez(k=size_z-2) вне домена — в маске");
    test::assert_eq(zero_violations, 0, "ZMAX: лишняя плоскость и Ez вне домена равны нулю");

    // Ez на k=size_z-3 лежит внутри домена — свободный DOF, не в маске.
    int ez_inside = 0;
    for (int i = 0; i < domain.size().x(); ++i) {
        for (int j = 0; j < domain.size().y(); ++j) {
            if (in_mask(domain, Face::ZMAX, i, j, k_face - 1, 2)) {
                ez_inside++;
            }
        }
    }
    test::assert_eq(ez_inside, 0, "ZMAX: Ez(k=size_z-3) не в маске");
}

void test_overwrite_has_no_mask() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    handler.load_from_json(
        bc_json(R"([{"er0": {"face": "ZMIN", "inner_radius": 1.0, "width": 1.0, "potential_drop": 0.1}}])"), domain);
    test::assert_true(handler.electric_dirichlet(domain).empty(), "overwrite: маска пуста по умолчанию");
}

void test_mask_respects_cylinder() {
    Domain domain = make_domain(true);
    BoundaryConditionHandler handler;
    handler.load_from_json(bc_json(R"([{"open": [{"face": "CYLINDER"}, {"face": "ZMAX"}]}, )"
                                   R"({"er0": {"face": "ZMIN", "inner_radius": 1.0, "width": 1.0, )"
                                   R"("potential_drop": 0.1, "mode": "solv"}}])"),
                           domain);
    const auto& entries = handler.electric_dirichlet(domain);
    int outside = 0;
    for (const auto& e : entries) {
        const auto [i, j, k, d] = domain.grid.pos_vind_range(e.row);
        const Vector3R pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
        if (std::hypot(pos.x() - 4.0, pos.y() - 4.0) > 4.0 + kMaskEps) {
            outside++;
        }
    }
    test::assert_true(!entries.empty(), "цилиндр: маска не пуста");
    test::assert_eq(outside, 0, "цилиндр: узлы вне цилиндра исключены");
}

void test_open_plus_er0_conflict() {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    bool thrown = false;
    try {
        handler.load_from_json(bc_json(R"([{"open": [{"face": "ZMIN"}]}, {"er0": {"face": "ZMIN", )"
                                       R"("inner_radius": 1.0, "width": 1.0, "potential_drop": 0.1, )"
                                       R"("mode": "solv"}}])"),
                               domain);
    } catch (const std::runtime_error&) {
        thrown = true;
    }
    test::assert_true(thrown, "open+er0 на одной грани -> ошибка");
}

bool load_throws(const std::string& conditions) {
    Domain domain = make_domain();
    BoundaryConditionHandler handler;
    try {
        handler.load_from_json(bc_json(conditions), domain);
    } catch (const std::exception&) {
        return true;
    }
    return false;
}

void test_invalid_er0_params() {
    test::assert_true(load_throws(R"([{"er0": {"face": "ZMIN", "width": 0.0}}])"), "er0: width=0 -> ошибка");
    test::assert_true(load_throws(R"([{"er0": {"face": "ZMIN", "inner_radius": -1.0}}])"),
                      "er0: inner_radius<0 -> ошибка");
    test::assert_true(load_throws(R"([{"er0": {"face": "XMIN"}}])"), "er0: face XMIN -> ошибка");
    test::assert_true(load_throws(R"([{"er0": {"face": "ZMIN", "mode": "fast"}}])"),
                      "er0: неизвестный mode -> ошибка");
    test::assert_true(!load_throws(R"([{"er0": {"face": "ZMIN"}}])"), "er0: параметры по умолчанию допустимы");
}

void test_impose_dirichlet_rows() {
    Operator mat(6, 6);
    std::vector<Trip> trips;
    for (int i = 0; i < 6; ++i) {
        trips.emplace_back(i, i, 1.0);
        if (i + 1 < 6) {
            trips.emplace_back(i, i + 1, 0.5);
            trips.emplace_back(i + 1, i, 0.5);
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
    mat.makeCompressed();

    impose_dirichlet_rows(mat, {{0, 2.5}, {3, -1.0}});

    test::assert_eq(mat.coeff(0, 0), 1.0, kEps, "impose_rows: A(0,0)=1");
    test::assert_eq(mat.coeff(0, 1), 0.0, kEps, "impose_rows: строка 0 обнулена");
    test::assert_eq(mat.coeff(3, 3), 1.0, kEps, "impose_rows: A(3,3)=1");
    test::assert_eq(mat.coeff(3, 2), 0.0, kEps, "impose_rows: строка 3 обнулена");
    // Столбцы закреплённых DOF не трогаются: свободные строки сохраняют связь
    // с предписанным полем (вклад A_fc * g).
    test::assert_eq(mat.coeff(1, 0), 0.5, kEps, "impose_rows: столбец 0 сохранён в свободной строке");
    test::assert_eq(mat.coeff(4, 3), 0.5, kEps, "impose_rows: столбец 3 сохранён в свободной строке");
    test::assert_eq(mat.coeff(1, 1), 1.0, kEps, "impose_rows: свободная строка не изменена");

    // Строка без диагонального элемента -> исключение.
    Operator bad(2, 2);
    std::vector<Trip> bad_trips = {Trip(0, 1, 1.0), Trip(1, 0, 1.0)};
    bad.setFromTriplets(bad_trips.begin(), bad_trips.end());
    bool thrown = false;
    try {
        impose_dirichlet_rows(bad, {{0, 1.0}});
    } catch (const std::runtime_error&) {
        thrown = true;
    }
    test::assert_true(thrown, "impose_rows: нет диагонали -> ошибка");
}

void test_impose_dirichlet_rhs() {
    Domain domain = make_domain();
    Field3d rhs(domain.size(), 3);
    rhs.setZero();

    impose_dirichlet_rhs(rhs, {{0, 2.5}, {3, -1.0}}, domain);

    const auto [i0, j0, k0, d0] = domain.grid.pos_vind_range(0);
    test::assert_eq(rhs(i0, j0, k0, d0), 2.5, kEps, "impose_rhs: строка 0");
    const auto [i1, j1, k1, d1] = domain.grid.pos_vind_range(3);
    test::assert_eq(rhs(i1, j1, k1, d1), -1.0, kEps, "impose_rhs: строка 3");

    double untouched = 0.0;
    for (int i = 0; i < domain.size().x(); ++i) {
        for (int j = 0; j < domain.size().y(); ++j) {
            for (int k = 0; k < domain.size().z(); ++k) {
                for (int d = 0; d < 3; ++d) {
                    const int row = domain.vind(i, j, k, d);
                    if (row != 0 && row != 3) {
                        untouched = std::max(untouched, std::abs(rhs(i, j, k, d)));
                    }
                }
            }
        }
    }
    test::assert_eq(untouched, 0.0, kEps, "impose_rhs: остальные значения не тронуты");
}

}   // namespace

int main() {
    test_mask_zmin();
    test_mask_zmax();
    test_overwrite_has_no_mask();
    test_mask_respects_cylinder();
    test_open_plus_er0_conflict();
    test_invalid_er0_params();
    test_impose_dirichlet_rows();
    test_impose_dirichlet_rhs();

    test::print_summary();

    return test::get_failed();
}
