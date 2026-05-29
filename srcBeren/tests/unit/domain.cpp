#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "vector3.h"
#include "World.h"

// Простой тестовый фреймворк без сторонних библиотек
namespace test {

static int passed = 0;
static int failed = 0;

void assert_true(bool condition, const std::string& test_name) {
    if (condition) {
      //  std::cout << "[PASS] " << test_name << std::endl;
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << std::endl;
        failed++;
    }
}

void assert_eq(double actual, double expected, double eps, const std::string& test_name) {
    if (std::abs(actual - expected) < eps) {
    //    std::cout << "[PASS] " << test_name << std::endl;
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << " (expected: " << expected 
                  << ", got: " << actual << ")" << std::endl;
        failed++;
    }
}

void assert_eq(int actual, int expected, const std::string& test_name) {
    if (actual == expected) {
     //   std::cout << "[PASS] " << test_name << std::endl;
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << " (expected: " << expected 
                  << ", got: " << actual << ")" << std::endl;
        failed++;
    }
}

void print_summary() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void reset() {
    passed = 0;
    failed = 0;
}

int get_failed() {
    return failed;
}

} // namespace test

// =============================================================================
// Тесты для Vector3R
// =============================================================================
void test_vector3_construction() {
    std::cout << "\n--- Vector3 Construction Tests ---" << std::endl;
    
    Vector3R v1;
    test::assert_eq(v1.x(), 0.0, 1e-10, "Vector3R default constructor x");
    test::assert_eq(v1.y(), 0.0, 1e-10, "Vector3R default constructor y");
    test::assert_eq(v1.z(), 0.0, 1e-10, "Vector3R default constructor z");
    
    Vector3R v2(1.0, 2.0, 3.0);
    test::assert_eq(v2.x(), 1.0, 1e-10, "Vector3R param constructor x");
    test::assert_eq(v2.y(), 2.0, 1e-10, "Vector3R param constructor y");
    test::assert_eq(v2.z(), 3.0, 1e-10, "Vector3R param constructor z");
    
    Vector3R v3(5.0);
    test::assert_eq(v3.x(), 5.0, 1e-10, "Vector3R uniform constructor x");
    test::assert_eq(v3.y(), 5.0, 1e-10, "Vector3R uniform constructor y");
    test::assert_eq(v3.z(), 5.0, 1e-10, "Vector3R uniform constructor z");
}

void test_vector3_operations() {
    std::cout << "\n--- Vector3 Operations Tests ---" << std::endl;
    
    Vector3R v1(1.0, 2.0, 3.0);
    Vector3R v2(4.0, 5.0, 6.0);
    
    Vector3R sum = v1 + v2;
    test::assert_eq(sum.x(), 5.0, 1e-10, "Vector3 addition x");
    test::assert_eq(sum.y(), 7.0, 1e-10, "Vector3 addition y");
    test::assert_eq(sum.z(), 9.0, 1e-10, "Vector3 addition z");
    
    Vector3R diff = v2 - v1;
    test::assert_eq(diff.x(), 3.0, 1e-10, "Vector3 subtraction x");
    test::assert_eq(diff.y(), 3.0, 1e-10, "Vector3 subtraction y");
    test::assert_eq(diff.z(), 3.0, 1e-10, "Vector3 subtraction z");
    
    Vector3R scaled = v1 * 2.0;
    test::assert_eq(scaled.x(), 2.0, 1e-10, "Vector3 scalar mult x");
    test::assert_eq(scaled.y(), 4.0, 1e-10, "Vector3 scalar mult y");
    test::assert_eq(scaled.z(), 6.0, 1e-10, "Vector3 scalar mult z");
    
    double dot = v1.dot(v2);
    test::assert_eq(dot, 32.0, 1e-10, "Vector3 dot product");
    
    Vector3R cross = v1.cross(v2);
    test::assert_eq(cross.x(), -3.0, 1e-10, "Vector3 cross product x");
    test::assert_eq(cross.y(), 6.0, 1e-10, "Vector3 cross product y");
    test::assert_eq(cross.z(), -3.0, 1e-10, "Vector3 cross product z");
}

// =============================================================================
// Тесты для Grid
// =============================================================================
void test_grid_construction() {
    std::cout << "\n--- Grid Construction Tests ---" << std::endl;
    
    Vector3R cell_size(0.5, 0.5, 0.5);
    Vector3I num_cells(10, 10, 10);
    int ghost_cells = 2;
    
    Grid grid(cell_size, num_cells, ghost_cells);
    
    test::assert_eq(grid.cell_size().x(), 0.5, 1e-10, "Grid cell_size x");
    test::assert_eq(grid.cell_size().y(), 0.5, 1e-10, "Grid cell_size y");
    test::assert_eq(grid.cell_size().z(), 0.5, 1e-10, "Grid cell_size z");
    
    test::assert_eq(grid.num_cells().x(), 10, "Grid num_cells x");
    test::assert_eq(grid.num_cells().y(), 10, "Grid num_cells y");
    test::assert_eq(grid.num_cells().z(), 10, "Grid num_cells z");
    
    test::assert_eq(grid.ghost_cells(), 2, "Grid ghost_cells");
    
    // size = num_cells + 2*ghost + 1 = 10 + 4 + 1 = 15
    test::assert_eq(grid.size().x(), 15, "Grid size x");
    test::assert_eq(grid.size().y(), 15, "Grid size y");
    test::assert_eq(grid.size().z(), 15, "Grid size z");
    
    double expected_volume = 0.5 * 0.5 * 0.5;
    test::assert_eq(grid.cell_volume(), expected_volume, 1e-10, "Grid cell_volume");
}

void test_grid_indexing() {
    std::cout << "\n--- Grid Indexing Tests ---" << std::endl;
    
    Vector3R cell_size(1.0, 1.0, 1.0);
    Vector3I num_cells(3, 4, 5); 
    int ghost_cells = 1;
    // size 6 7 8
    Grid grid(cell_size, num_cells, ghost_cells);
    
    // sind(i, j, k) = i * size.y * size.z + j * size.z + k
    // 2, 3, 1
    int expected_sind = 2 * 7 * 8 + 3 * 8 + 1;
    test::assert_eq(grid.sind(2, 3, 1), expected_sind, "Grid sind");
    
    // vind(i, j, k, d, nd) = d + nd * (i * size.y * size.z + j * size.z + k)
    int expected_vind = 2 + 3 * (1 * 7 * 8 + 2 * 8 + 3);
    test::assert_eq(grid.vind(1, 2, 3, 2), expected_vind, "Grid vind");
}

void test_grid_coordinates() {
    std::cout << "\n--- Grid Coordinate Tests ---" << std::endl;
    
    Vector3R cell_size(0.5, 0.5, 0.5);
    Vector3I num_cells(11, 12, 13);
    int ghost_cells = 1;
    
    Grid grid(cell_size, num_cells, ghost_cells);
    
    Vector3R world(1.0, 2.0, 3.0);
    Vector3R cell = grid.to_cell_coordinates(world);
    test::assert_eq(cell.x(), 2.0, 1e-10, "Grid to_cell_coordinates x");
    test::assert_eq(cell.y(), 4.0, 1e-10, "Grid to_cell_coordinates y");
    test::assert_eq(cell.z(), 6.0, 1e-10, "Grid to_cell_coordinates z");
    
    Vector3R coord(1.25, 1.75, 2.0);
    Vector3I idx = grid.get_cell_index(coord);
    // idx = int(coord / cell_size + ghost) = int(2.5 + 1) = 3
    // idy = int(coord / cell_size + ghost) = int(3.5 + 1) = 4
    // idz = int(coord / cell_size + ghost) = int(4 + 1) = 5
    test::assert_eq(idx.x(), 3, "Grid get_cell_index x");
    test::assert_eq(idx.y(), 4, "Grid get_cell_index y");
    test::assert_eq(idx.z(), 5, "Grid get_cell_index z");
}

// =============================================================================
// Тесты для Geometry
// =============================================================================
void test_geometry_box() {
    std::cout << "\n--- Geometry Box Tests ---" << std::endl;
    
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(11.0, 13.0, 12.0);
    
    Geometry geom(box_min, box_max, false);
    
    Vector3R inside(5.0, 5.0, 5.0);
    test::assert_true(geom.contains(inside, 1e-10), "Geometry contains point inside box");
    
    Vector3R outside(-0.01, 3.0, 12.0);
    test::assert_true(!geom.contains(outside, 1e-10), "Geometry rejects point outside box (x < min)");
    
    Vector3R outside2(11.0, 5.0, 5.0);
    test::assert_true(!geom.contains(outside2), "Geometry rejects point on boundary (x = max)");
    
    Vector3R on_min(0.001, 0.0001, 0.00001);
    test::assert_true(geom.contains(on_min, 1e-10), "Geometry contains point near min boundary");
}

void test_geometry_cylinder() {
    std::cout << "\n--- Geometry Cylinder Tests ---" << std::endl;
    
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(13.0, 13.0, 12.0);
    
    Geometry geom(box_min, box_max, true);
    geom.cyl_center = Vector3R(6.5, 6.5, 0.0);
    geom.cyl_radius = 6.5;
    
    // Точка внутри цилиндра (на оси)
    Vector3R inside_cyl(5.0, 5.0, 5.0);
    test::assert_true(geom.contains(inside_cyl, 1e-10), "Geometry contains point inside cylinder");
    
    // Точка вне цилиндра (слишком далеко от оси)
    Vector3R outside_cyl(-1.5, 1.5, 5.0);
    test::assert_true(!geom.contains(outside_cyl, 1e-10), "Geometry rejects point outside cylinder");
    
    // Точка на границе цилиндра
    Vector3R on_boundary(0.0, 0.0, 1.0);
    test::assert_true(!geom.contains(on_boundary, 1e-10), "Geometry rejects point on cylinder boundary");
    
    // Точка внутри box, но вне цилиндра по x
    Vector3R outside_box(1.5, 1.5, 5.0);
    // distance^2 = (1.5-6.5)^2 + (1.5-6.5)^2 = 25 + 25 = 50 > 6.5*6.5
    test::assert_true(!geom.contains(outside_box, 1e-10), "Geometry rejects point outside box");
}

void test_geometry_faces() {
    std::cout << "\n--- Geometry Face Tests ---" << std::endl;
    
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(8.0, 13.0, 20.0);
    
    Geometry geom(box_min, box_max, false);
    
    Vector3R p(8.0, 3.0, 6.0);
    test::assert_true(!geom.is_outside_face(Face::XMIN, p, 1e-10), "Point not outside XMIN");
    
    Vector3R p_xmin(-0.01, 3.0, 15.0);
    test::assert_true(geom.is_outside_face(Face::XMIN, p_xmin, 1e-10), "Point outside XMIN");
    
    Vector3R p_xmax(8.0, 2.0, 0.01);
    test::assert_true(geom.is_outside_face(Face::XMAX, p_xmax, 1e-10), "Point outside XMAX");
    
    Vector3R p_ymin(5.0, -0.00004, 5.0);
    test::assert_true(geom.is_outside_face(Face::YMIN, p_ymin, 1e-10), "Point outside YMIN");
    
    Vector3R p_ymax(-1.0, 13.0, 5.0);
    test::assert_true(geom.is_outside_face(Face::YMAX, p_ymax, 1e-10), "Point outside YMAX");
    
    Vector3R p_zmin(-5.0, 5.0, -0.0001);
    test::assert_true(geom.is_outside_face(Face::ZMIN, p_zmin, 1e-10), "Point outside ZMIN");
    
    Vector3R p_zmax(5.0, 5.0, 20.0000000001);
    test::assert_true(geom.is_outside_face(Face::ZMAX, p_zmax, 1e-10), "Point outside ZMAX");
}

// =============================================================================
// Тесты для Domain
// =============================================================================
void test_domain_init() {
    std::cout << "\n--- Domain Initialization Tests ---" << std::endl;
    
    Domain domain;
    Vector3I num_cells(105, 10, 5);
    Vector3R cell_size(0.1, 0.2, 0.3);
    int ghost = 3;
    
    domain.init(num_cells, cell_size, ghost);
    
    test::assert_eq(domain.grid.cell_size().x(), 0.1, 1e-10, "Domain grid cell_size x");
    test::assert_eq(domain.grid.num_cells().x(), 105, "Domain grid num_cells x");
    test::assert_eq(domain.grid.ghost_cells(), 3, "Domain grid ghost_cells");
    
    // box_max = num_cells * cell_size
    test::assert_eq(domain.geom.box_max.x(), 10.5, 1e-10, "Domain geom box_max x");
    test::assert_eq(domain.geom.box_max.y(), 2.0, 1e-10, "Domain geom box_max y");
    test::assert_eq(domain.geom.box_max.z(), 1.5, 1e-10, "Domain geom box_max z");
}

void test_domain_cylinder() {
    std::cout << "\n--- Domain Cylinder Tests ---" << std::endl;
    
    Domain domain;
    Vector3I num_cells(20, 20, 20);
    Vector3R cell_size(0.2, 0.2, 0.2);
    
    domain.init(num_cells, cell_size, 1);
    domain.set_cylinder(Vector3R(2.0, 2.0, 2.0), 2.0);
    
    test::assert_true(domain.geom.use_cylinder, "Domain cylinder enabled");
    test::assert_eq(domain.geom.cyl_center.x(), 2.0, 1e-10, "Domain cyl_center x");
    test::assert_eq(domain.geom.cyl_center.y(), 2.0, 1e-10, "Domain cyl_center y");
    test::assert_eq(domain.geom.cyl_radius, 2.0, 1e-10, "Domain cyl_radius");
    
    // Точка внутри цилиндра
    Vector3R inside(2.5, 1.5, 1.5);
    test::assert_true(domain.contains(inside, 1e-10), "Domain contains point inside cylinder");
    
    // Точка вне цилиндра
    Vector3R outside(3.8, 3.5, 3.0);
    test::assert_true(!domain.contains(outside, 1e-10), "Domain rejects point outside cylinder");
}

void test_domain_node_position() {
    std::cout << "\n--- Domain Node Position Tests ---" << std::endl;
    
    Domain domain;
    Vector3I num_cells(12, 13, 11);
    Vector3R cell_size(0.5, 0.5, 0.5);
    
    domain.init(num_cells, cell_size, 1);
    
    // Тест для электрического поля (Ex компонент)
    // Ex сдвинута на +dx/2 по X
    Vector3R pos_ex = domain.get_node_position(1, 1, 1, FieldType::ELECTRIC, 0);
    // x0 = origin.x + (i - ghost) * cell_size = 0 + (1 - 1) * 1 = 0
    // pos_x = x0 + 0.5 * cell_size = 0 + 0.5 = 0.5
    test::assert_eq(pos_ex.x(), 0.25, 1e-10, "Domain Ex node position x");
    test::assert_eq(pos_ex.y(), 0.0, 1e-10, "Domain Ex node position y");
    test::assert_eq(pos_ex.z(), 0.0, 1e-10, "Domain Ex node position z");

    Vector3R pos_ey = domain.get_node_position(0, 2, 1, FieldType::ELECTRIC, 1);
    test::assert_eq(pos_ey.x(), -0.5, 1e-10, "Domain Ex node position x");
    test::assert_eq(pos_ey.y(), 0.75, 1e-10, "Domain Ex node position y");
    test::assert_eq(pos_ey.z(), 0.0, 1e-10, "Domain Ex node position z");

    Vector3R pos_ez = domain.get_node_position(0, 1, 1, FieldType::ELECTRIC, 2);
    test::assert_eq(pos_ez.x(), -0.5, 1e-10, "Domain Ex node position x");
    test::assert_eq(pos_ez.y(), 0.0, 1e-10, "Domain Ex node position y");
    test::assert_eq(pos_ez.z(), 0.25, 1e-10, "Domain Ex node position z");

    // Тест для магнитного поля (Bx компонент)
    // Bx сдвинута на +dy/2 и +dz/2
    Vector3R pos_bx = domain.get_node_position(1, 1, 1, FieldType::MAGNETIC, 0);
    test::assert_eq(pos_bx.x(), 0.0, 1e-10, "Domain Bx node position x");
    test::assert_eq(pos_bx.y(), 0.25, 1e-10, "Domain Bx node position y");
    test::assert_eq(pos_bx.z(), 0.25, 1e-10, "Domain Bx node position z");
    Vector3R pos_by = domain.get_node_position(1, 1, 1, FieldType::MAGNETIC, 1);
    test::assert_eq(pos_by.x(), 0.25, 1e-10, "Domain Bx node position x");
    test::assert_eq(pos_by.y(), 0.0, 1e-10, "Domain Bx node position y");
    test::assert_eq(pos_by.z(), 0.25, 1e-10, "Domain Bx node position z");
    Vector3R pos_bz = domain.get_node_position(1, 1, 1, FieldType::MAGNETIC, 2);
    test::assert_eq(pos_bz.x(), 0.25, 1e-10, "Domain Bx node position x");
    test::assert_eq(pos_bz.y(), 0.25, 1e-10, "Domain Bx node position y");
    test::assert_eq(pos_bz.z(), 0.0, 1e-10, "Domain Bx node position z");
    // Тест для плотности (в центре ячейки)
    Vector3R pos_density = domain.get_node_position(1, 1, 1, FieldType::DENSITY, 0);
    test::assert_eq(pos_density.x(), 0.0, 1e-10, "Domain density node position x");
    test::assert_eq(pos_density.y(), 0.0, 1e-10, "Domain density node position y");
    test::assert_eq(pos_density.z(), 0.0, 1e-10, "Domain density node position z");
}

void test_domain_is_inside_node() {
    std::cout << "\n--- Domain Is Inside Node Tests ---" << std::endl;
    
    Domain domain;
    Vector3I num_cells(14, 14, 7);
    Vector3R cell_size(0.5, 0.5, 0.5);
    
    domain.init(num_cells, cell_size, 1);
    auto size = domain.grid.size();

    for (int i = 0; i < size.x(); i++){
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                bool inside_ex =
                    domain.is_inside_node(i, j, k, FieldType::ELECTRIC, 0);
                    if( i < 1 || j < 2 || k < 2 || i > size.x() - 3 || j > size.y() -3 || k > size.z() -3){
                        test::assert_true(!inside_ex,
                                          "Domain node outside box");
                    } else {
                        test::assert_true(inside_ex, "Domain node inside box");
                    }
                    bool inside_ey =
                        domain.is_inside_node(i, j, k, FieldType::ELECTRIC, 1);
                    if (i < 2 || j < 1 || k < 2 || i > size.x() - 3 ||
                        j > size.y() - 3 || k > size.z() - 3) {
                        test::assert_true(!inside_ey,
                                          "Domain node outside box");
                    } else {
                        test::assert_true(inside_ey, "Domain node inside box");
                    }
                    bool inside_ez =
                        domain.is_inside_node(i, j, k, FieldType::ELECTRIC, 2);
                    if (i < 2 || j < 2 || k < 1 || i > size.x() - 3 ||
                        j > size.y() - 3 || k > size.z() - 3) {
                        test::assert_true(!inside_ez,
                                          "Domain node outside box");
                    } else {
                        test::assert_true(inside_ez, "Domain node inside box");
                    }
                    bool inside_bx =
                        domain.is_inside_node(i, j, k, FieldType::MAGNETIC, 0);
                    if (i < 2 || j < 1 || k < 1 || i > size.x() - 3 ||
                        j > size.y() - 3 || k > size.z() - 3) {
                        test::assert_true(!inside_bx,
                                          "Domain node outside box");
                    } else {
                        test::assert_true(inside_bx, "Domain node inside box");
                    }
                    bool inside_by =
                        domain.is_inside_node(i, j, k, FieldType::MAGNETIC, 1);
                    if (i < 1 || j < 2 || k < 1 || i > size.x() - 3 ||
                        j > size.y() - 3 || k > size.z() - 3) {
                        test::assert_true(!inside_by,
                                          "Domain node outside box");
                    } else {
                        test::assert_true(inside_by, "Domain node inside box");
                    }
                    bool inside_bz =
                        domain.is_inside_node(i, j, k, FieldType::MAGNETIC, 2);
                    if (i < 1 || j < 1 || k < 2 || i > size.x() - 3 ||
                        j > size.y() - 3 || k > size.z() - 3) {
                        test::assert_true(!inside_bz,
                                          "Domain node outside box");
                    } else {
                        test::assert_true(inside_bz, "Domain node inside box");
                    }
            }
        }
    }

    domain.set_cylinder(Vector3R(3.5, 3.5, 0.0), 3.5);
    
    // Точка внутри цилиндра
    bool inside = domain.is_inside_node(5, 6, 2, FieldType::DENSITY, 0);
    test::assert_true(inside, "Domain node inside cylinder");
    
    // Точка вне цилиндра
    bool outside = domain.is_inside_node(2, 3, 3, FieldType::DENSITY, 0);
    test::assert_true(!outside, "Domain node outside cylinder");
}

// =============================================================================
// Дополнительные тесты для Grid
// =============================================================================
void test_grid_origin_and_total() {
    std::cout << "\n--- Grid Origin & Total Size Tests ---" << std::endl;

    Vector3R cell_size(1.0, 2.0, 3.0);
    Vector3I num_cells(5, 6, 7);
    int ghost_cells = 2;

    Grid grid(cell_size, num_cells, ghost_cells);

    test::assert_eq(grid.origin().x(), 0.0, 1e-10, "Grid origin x");
    test::assert_eq(grid.origin().y(), 0.0, 1e-10, "Grid origin y");
    test::assert_eq(grid.origin().z(), 0.0, 1e-10, "Grid origin z");

    // size = (5+4+1, 6+4+1, 7+4+1) = (10, 11, 12)
    int expected_total = 10 * 11 * 12;
    test::assert_eq(grid.total_size(), expected_total, "Grid total_size");
}

void test_grid_pos_vind() {
    std::cout << "\n--- Grid pos_vind Tests ---" << std::endl;

    Vector3R cell_size(1.0, 1.0, 1.0);
    Vector3I num_cells(3, 4, 5);
    int ghost_cells = 1;
    // size: 6, 7, 8; dims_ = {6, 7, 8, 3}

    Grid grid(cell_size, num_cells, ghost_cells);

    int index = 0;
    test::assert_eq(grid.pos_vind(index, 0), 0, "Grid pos_vind index=0 n=0");
    test::assert_eq(grid.pos_vind(index, 1), 0, "Grid pos_vind index=0 n=1");
    test::assert_eq(grid.pos_vind(index, 2), 0, "Grid pos_vind index=0 n=2");
    test::assert_eq(grid.pos_vind(index, 3), 0, "Grid pos_vind index=0 n=3");

    // index = d + nd*(i*sy*sz + j*sz + k) = 2 + 3*(1*7*8 + 2*8 + 3) = 2 + 3*(56+16+3) = 2 + 225 = 227
    index = 227;
    test::assert_eq(grid.pos_vind(index, 0), 1, "Grid pos_vind index=227 n=0 (i)");
    test::assert_eq(grid.pos_vind(index, 1), 2, "Grid pos_vind index=227 n=1 (j)");
    test::assert_eq(grid.pos_vind(index, 2), 3, "Grid pos_vind index=227 n=2 (k)");
    test::assert_eq(grid.pos_vind(index, 3), 2, "Grid pos_vind index=227 n=3 (d)");
}

void test_grid_pos_sind() {
    std::cout << "\n--- Grid pos_sind Tests ---" << std::endl;

    Vector3R cell_size(1.0, 1.0, 1.0);
    Vector3I num_cells(3, 4, 5);
    int ghost_cells = 1;
    // size: 6, 7, 8

    Grid grid(cell_size, num_cells, ghost_cells);

    int index = 0;
    test::assert_eq(grid.pos_sind(index, 0), 0, "Grid pos_sind index=0 n=0");
    test::assert_eq(grid.pos_sind(index, 1), 0, "Grid pos_sind index=0 n=1");
    test::assert_eq(grid.pos_sind(index, 2), 0, "Grid pos_sind index=0 n=2");

    // index = i*sy*sz + j*sz + k = 1*7*8 + 2*8 + 3 = 56 + 16 + 3 = 75
    index = 75;
    test::assert_eq(grid.pos_sind(index, 0), 1, "Grid pos_sind index=75 n=0 (i)");
    test::assert_eq(grid.pos_sind(index, 1), 2, "Grid pos_sind index=75 n=1 (j)");
    test::assert_eq(grid.pos_sind(index, 2), 3, "Grid pos_sind index=75 n=2 (k)");
}

void test_grid_get_field_node_index() {
    std::cout << "\n--- Grid get_field_node_index Tests ---" << std::endl;

    Vector3R cell_size(0.5, 0.5, 0.5);
    Vector3I num_cells(10, 10, 10);
    int ghost_cells = 1;

    Grid grid(cell_size, num_cells, ghost_cells);
    Vector3R world(1.0, 2.0, 3.0);

    // Ex: shift = (0.5, 0, 0) in cell units
    // x = 1.0/0.5 - 0.5 = 2.0 - 0.5 = 1.5, floor = 1, i = 1 + 1 = 2
    // y = 2.0/0.5 - 0 = 4.0, floor = 4, j = 4 + 1 = 5
    // z = 3.0/0.5 - 0 = 6.0, floor = 6, k = 6 + 1 = 7
    Vector3I idx_ex = grid.get_field_node_index(world, FieldType::ELECTRIC, 0);
    test::assert_eq(idx_ex.x(), 2, "Grid Ex node index x");
    test::assert_eq(idx_ex.y(), 5, "Grid Ex node index y");
    test::assert_eq(idx_ex.z(), 7, "Grid Ex node index z");

    // Bx: shift = (0, 0.5, 0.5) in cell units
    // x = 2.0 - 0 = 2.0, floor = 2, i = 3
    // y = 4.0 - 0.5 = 3.5, floor = 3, j = 4
    // z = 6.0 - 0.5 = 5.5, floor = 5, k = 6
    Vector3I idx_bx = grid.get_field_node_index(world, FieldType::MAGNETIC, 0);
    test::assert_eq(idx_bx.x(), 3, "Grid Bx node index x");
    test::assert_eq(idx_bx.y(), 4, "Grid Bx node index y");
    test::assert_eq(idx_bx.z(), 6, "Grid Bx node index z");

    // Density: shift = (0, 0, 0)
    // x = 2.0, floor = 2, i = 3
    // y = 4.0, floor = 4, j = 5
    // z = 6.0, floor = 6, k = 7
    Vector3I idx_rho = grid.get_field_node_index(world, FieldType::DENSITY, 0);
    test::assert_eq(idx_rho.x(), 3, "Grid density node index x");
    test::assert_eq(idx_rho.y(), 5, "Grid density node index y");
    test::assert_eq(idx_rho.z(), 7, "Grid density node index z");
}

// =============================================================================
// Дополнительные тесты для Geometry
// =============================================================================
void test_geometry_in_cylinder() {
    std::cout << "\n--- Geometry in_cylinder Tests ---" << std::endl;

    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);

    Geometry geom;
    geom.init(box_min, box_max, false);
    test::assert_true(!geom.in_cylinder(Vector3R(5.0, 5.0, 5.0), 1e-10),
                      "in_cylinder false when use_cylinder=false");

    geom.use_cylinder = true;
    geom.cyl_center = Vector3R(5.0, 5.0, 0.0);
    geom.cyl_radius = 3.0;

    Vector3R inside(5.0, 5.0, 5.0);
    test::assert_true(geom.in_cylinder(inside, 1e-10), "in_cylinder at center");

    Vector3R near_boundary(5.0, 7.9, 5.0);
    test::assert_true(geom.in_cylinder(near_boundary, 1e-10), "in_cylinder near boundary");

    Vector3R outside(5.0, 8.1, 5.0);
    test::assert_true(!geom.in_cylinder(outside, 1e-10), "in_cylinder outside");
}

void test_geometry_contains_periodic() {
    std::cout << "\n--- Geometry contains(periodic) Tests ---" << std::endl;

    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    bool periodic[3] = {false, false, false};
    // x outside
    test::assert_true(!geom.contains(Vector3R(-1.0, 5.0, 5.0), periodic, 1e-10),
                      "contains periodic: x < min rejected");
    // periodic x — пропускает границу
    bool periodic_x[3] = {true, false, false};
    test::assert_true(geom.contains(Vector3R(-1.0, 5.0, 5.0), periodic_x, 1e-10),
                      "contains periodic: x periodic skips x min");
    // periodic all — любая точка внутри
    bool periodic_all[3] = {true, true, true};
    test::assert_true(geom.contains(Vector3R(-999.0, 999.0, -999.0), periodic_all, 1e-10),
                      "contains periodic: all periodic accepts anything");
}

void test_geometry_contains_ignoring_face() {
    std::cout << "\n--- Geometry contains_ignoring_face Tests ---" << std::endl;

    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Точка вне XMIN, но XMIN игнорируется — внутри
    Vector3R p(-1.0, 5.0, 5.0);
    test::assert_true(geom.contains_ignoring_face(Face::XMIN, p, 1e-10),
                      "contains_ignoring_face XMIN");
    test::assert_true(!geom.contains_ignoring_face(Face::XMAX, p, 1e-10),
                      "contains_ignoring_face XMAX rejects x < min");

    // Точка вне XMAX, но XMAX игнорируется — внутри
    Vector3R p2(10.5, 5.0, 5.0);
    test::assert_true(geom.contains_ignoring_face(Face::XMAX, p2, 1e-10),
                      "contains_ignoring_face XMAX");

    // CYLINDER
    geom.use_cylinder = true;
    geom.cyl_center = Vector3R(5.0, 5.0, 0.0);
    geom.cyl_radius = 3.0;
    Vector3R outside_cyl(5.0, 9.0, 5.0);
    test::assert_true(!geom.contains_ignoring_face(Face::XMIN, outside_cyl, 1e-10),
                      "contains_ignoring_face rejects outside cylinder");
    test::assert_true(geom.contains_ignoring_face(Face::CYLINDER, outside_cyl, 1e-10),
                      "contains_ignoring_face CYLINDER ignores cylinder check");
}

void test_geometry_is_outside_only_face() {
    std::cout << "\n--- Geometry is_outside_only_face Tests ---" << std::endl;

    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Точка чуть левее XMIN, все остальные грани внутри
    Vector3R p(-0.1, 5.0, 5.0);
    test::assert_true(geom.is_outside_only_face(Face::XMIN, p, 1e-10),
                      "is_outside_only_face XMIN only");

    // Точка вне XMIN и YMIN одновременно
    Vector3R p2(-0.1, -0.1, 5.0);
    test::assert_true(!geom.is_outside_only_face(Face::XMIN, p2, 1e-10),
                      "is_outside_only_face: outside two faces, not only XMIN");
}

void test_geometry_reflect_from_face() {
    std::cout << "\n--- Geometry reflect_from_face Tests ---" << std::endl;

    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Отражаем от XMIN: x' = 2*0 - (-1) = 1
    Vector3R p(-1.0, 5.0, 5.0);
    Vector3R refl = geom.reflect_from_face(Face::XMIN, p);
    test::assert_eq(refl.x(), 1.0, 1e-10, "reflect_from_face XMIN x");
    test::assert_eq(refl.y(), 5.0, 1e-10, "reflect_from_face XMIN y");
    test::assert_eq(refl.z(), 5.0, 1e-10, "reflect_from_face XMIN z");

    // Отражаем от XMAX: x' = 2*10 - 11 = 9
    Vector3R p2(11.0, 3.0, 7.0);
    Vector3R refl2 = geom.reflect_from_face(Face::XMAX, p2);
    test::assert_eq(refl2.x(), 9.0, 1e-10, "reflect_from_face XMAX x");

    // YMIN
    Vector3R p3(5.0, -2.0, 5.0);
    Vector3R refl3 = geom.reflect_from_face(Face::YMIN, p3);
    test::assert_eq(refl3.y(), 2.0, 1e-10, "reflect_from_face YMIN y");

    // YMAX
    Vector3R p4(5.0, 12.0, 5.0);
    Vector3R refl4 = geom.reflect_from_face(Face::YMAX, p4);
    test::assert_eq(refl4.y(), 8.0, 1e-10, "reflect_from_face YMAX y");

    // ZMIN
    Vector3R p5(5.0, 5.0, -3.0);
    Vector3R refl5 = geom.reflect_from_face(Face::ZMIN, p5);
    test::assert_eq(refl5.z(), 3.0, 1e-10, "reflect_from_face ZMIN z");

    // ZMAX
    Vector3R p6(5.0, 5.0, 13.0);
    Vector3R refl6 = geom.reflect_from_face(Face::ZMAX, p6);
    test::assert_eq(refl6.z(), 7.0, 1e-10, "reflect_from_face ZMAX z");
}

void test_geometry_reflect_cylinder() {
    std::cout << "\n--- Geometry reflect CYLINDER Tests ---" << std::endl;

    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);
    geom.use_cylinder = true;
    geom.cyl_center = Vector3R(5.0, 5.0, 0.0);
    geom.cyl_radius = 4.0;

    // Точка на расстоянии radius+1 = 5 от центра по оси x
    // dx=5, new_dist = 2*4 - 5 = 3, scale = 3/5 = 0.6
    // reflected.x = 5 + 5*0.6 = 8
    Vector3R p(10.0, 5.0, 3.0);
    Vector3R refl = geom.reflect_from_face(Face::CYLINDER, p);
    test::assert_eq(refl.x(), 8.0, 1e-10, "reflect CYLINDER x");
    test::assert_eq(refl.y(), 5.0, 1e-10, "reflect CYLINDER y");
    test::assert_eq(refl.z(), 3.0, 1e-10, "reflect CYLINDER z (unchanged)");
}

void test_geometry_reflect_from_boundary() {
    std::cout << "\n--- Geometry reflect_from_boundary Tests ---" << std::endl;

    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Точка за XMIN
    Vector3R p(-1.0, 5.0, 5.0);
    Vector3R refl = geom.reflect_from_boundary(p);
    test::assert_eq(refl.x(), 1.0, 1e-10, "reflect_boundary xmin x");
    test::assert_true(geom.contains(refl, 1e-12), "reflect_boundary: result inside box");

    // Точка за двумя гранями (XMIN и YMAX)
    Vector3R p2(-2.0, 13.0, 5.0);
    Vector3R refl2 = geom.reflect_from_boundary(p2);
    test::assert_true(geom.contains(refl2, 1e-12), "reflect_boundary: result inside after double reflect");

    // Точка за XMAX, ZMIN
    Vector3R p3(12.0, 5.0, -1.0);
    Vector3R refl3 = geom.reflect_from_boundary(p3);
    test::assert_true(geom.contains(refl3, 1e-12), "reflect_boundary: triple axis check inside");
}

// =============================================================================
// Дополнительные тесты для Domain
// =============================================================================
void test_domain_accessors() {
    std::cout << "\n--- Domain Accessor Tests ---" << std::endl;

    Domain domain;
    Vector3I num_cells(4, 5, 6);
    Vector3R cell_size(0.5, 0.5, 0.5);

    domain.init(num_cells, cell_size, 1);
    // size: 7, 8, 9; total = 7*8*9 = 504

    test::assert_eq(domain.total_size(), 504, "Domain total_size");
    test::assert_eq(domain.cell_size(0), 0.5, 1e-10, "Domain cell_size(int) x");
    test::assert_eq(domain.cell_size(1), 0.5, 1e-10, "Domain cell_size(int) y");
    test::assert_eq(domain.cell_size(2), 0.5, 1e-10, "Domain cell_size(int) z");
    test::assert_eq(domain.num_cells(0), 4, "Domain num_cells(int) x");
    test::assert_eq(domain.num_cells(1), 5, "Domain num_cells(int) y");
    test::assert_eq(domain.num_cells(2), 6, "Domain num_cells(int) z");
    test::assert_eq(domain.cell_volume(), 0.125, 1e-10, "Domain cell_volume");
}

void test_domain_sind_vind() {
    std::cout << "\n--- Domain sind/vind Tests ---" << std::endl;

    Domain domain;
    Vector3I num_cells(3, 4, 5);
    Vector3R cell_size(1.0, 1.0, 1.0);

    domain.init(num_cells, cell_size, 1);
    // size: 6, 7, 8

    // sind: i*sy*sz + j*sz + k = 2*7*8 + 3*8 + 4 = 112 + 24 + 4 = 140
    test::assert_eq(domain.sind(2, 3, 4), 140, "Domain sind");

    // vind: d + nd*(i*sy*sz + j*sz + k) = 1 + 3*(1*7*8 + 2*8 + 3) = 1 + 3*(56+16+3) = 1 + 225 = 226
    test::assert_eq(domain.vind(1, 2, 3, 1), 226, "Domain vind");
}

void test_domain_get_cell_index() {
    std::cout << "\n--- Domain get_cell_index Tests ---" << std::endl;

    Domain domain;
    Vector3I num_cells(10, 10, 10);
    Vector3R cell_size(0.5, 0.5, 0.5);

    domain.init(num_cells, cell_size, 1);

    Vector3R coord(0.75, 1.25, 2.0);
    Vector3I idx = domain.get_cell_index(coord);
    // idx = int(0.75/0.5 + 1) = int(1.5+1) = 2
    test::assert_eq(idx.x(), 2, "Domain get_cell_index x");
    test::assert_eq(idx.y(), 3, "Domain get_cell_index y");
    test::assert_eq(idx.z(), 5, "Domain get_cell_index z");
}

void test_domain_to_cell_coordinates() {
    std::cout << "\n--- Domain to_cell_coordinates Tests ---" << std::endl;

    Domain domain;
    Vector3I num_cells(10, 10, 10);
    Vector3R cell_size(0.5, 0.5, 0.5);

    domain.init(num_cells, cell_size, 1);

    Vector3R world(1.0, 2.0, 3.0);
    Vector3R cell = domain.to_cell_coordinates(world);
    test::assert_eq(cell.x(), 2.0, 1e-10, "Domain to_cell_coordinates x");
    test::assert_eq(cell.y(), 4.0, 1e-10, "Domain to_cell_coordinates y");
    test::assert_eq(cell.z(), 6.0, 1e-10, "Domain to_cell_coordinates z");
}

void test_domain_pos_vind_pos_sind() {
    std::cout << "\n--- Domain pos_vind/pos_sind Tests ---" << std::endl;

    Domain domain;
    Vector3I num_cells(3, 4, 5);
    Vector3R cell_size(1.0, 1.0, 1.0);

    domain.init(num_cells, cell_size, 1);

    // index = 227 = d + nd*(i*sy*sz + j*sz + k) with nd=3, i=1, j=2, k=3, d=2
    int idx = 227;
    test::assert_eq(domain.pos_vind(idx, 0), 1, "Domain pos_vind n=0");
    test::assert_eq(domain.pos_vind(idx, 1), 2, "Domain pos_vind n=1");
    test::assert_eq(domain.pos_vind(idx, 2), 3, "Domain pos_vind n=2");
    test::assert_eq(domain.pos_vind(idx, 3), 2, "Domain pos_vind n=3");

    int sidx = 75;  // 1*7*8 + 2*8 + 3 = 56+16+3
    test::assert_eq(domain.pos_sind(sidx, 0), 1, "Domain pos_sind n=0");
    test::assert_eq(domain.pos_sind(sidx, 1), 2, "Domain pos_sind n=1");
    test::assert_eq(domain.pos_sind(sidx, 2), 3, "Domain pos_sind n=2");
}

void test_domain_is_inside_node_linear() {
    std::cout << "\n--- Domain is_inside_node(linear_index) Tests ---" << std::endl;

    Domain domain;
    Vector3I num_cells(5, 5, 5);
    Vector3R cell_size(1.0, 1.0, 1.0);

    domain.init(num_cells, cell_size, 1);

    // vind(3, 3, 3, 0) для ELECTRIC (Ex) — внутри бокса
    int idx_ex = domain.vind(3, 3, 3, 0);
    test::assert_true(domain.is_inside_node(idx_ex, FieldType::ELECTRIC),
                      "Domain is_inside_node linear Ex inside");
    // vind(0, 3, 3, 0) — снаружи
    int idx_out = domain.vind(0, 3, 3, 0);
    test::assert_true(!domain.is_inside_node(idx_out, FieldType::ELECTRIC),
                      "Domain is_inside_node linear Ex outside");
    // vind(3, 3, 3, 1) для ELECTRIC (Ey) — внутри
    int idx_ey = domain.vind(3, 3, 3, 1);
    test::assert_true(domain.is_inside_node(idx_ey, FieldType::ELECTRIC),
                      "Domain is_inside_node linear Ey inside");
}

void test_domain_is_inside_node_periodic() {
    std::cout << "\n--- Domain is_inside_node_periodic Tests ---" << std::endl;

    Domain domain;
    Vector3I num_cells(5, 5, 5);
    Vector3R cell_size(1.0, 1.0, 1.0);

    domain.init(num_cells, cell_size, 1);

    bool no_periodic[3] = {false, false, false};
    bool periodic_x[3] = {true, false, false};

    // Крайний узел Ex (i=0) — за левой границей без periodicity
    test::assert_true(!domain.is_inside_node_periodic(0, 3, 3, FieldType::ELECTRIC, 0, no_periodic),
                      "is_inside_node_periodic: Ex at i=0 non-periodic rejected");
    // С periodicity по X — принимается
    test::assert_true(domain.is_inside_node_periodic(0, 3, 3, FieldType::ELECTRIC, 0, periodic_x),
                      "is_inside_node_periodic: Ex at i=0 periodic x accepted");
}

// =============================================================================
// Главная функция
// =============================================================================
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Running Boundary Tests" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Vector3 тесты
    test_vector3_construction();
    test_vector3_operations();

    // Grid тесты
    test_grid_construction();
    test_grid_indexing();
    test_grid_coordinates();
    test_grid_origin_and_total();
    test_grid_pos_vind();
    test_grid_pos_sind();
    test_grid_get_field_node_index();

    // Geometry тесты
    test_geometry_box();
    test_geometry_cylinder();
    test_geometry_faces();
    test_geometry_in_cylinder();
    test_geometry_contains_periodic();
    test_geometry_contains_ignoring_face();
    test_geometry_is_outside_only_face();
    test_geometry_reflect_from_face();
    test_geometry_reflect_cylinder();
    test_geometry_reflect_from_boundary();

    // Domain тесты
    test_domain_init();
    test_domain_cylinder();
    test_domain_node_position();
    test_domain_is_inside_node();
    test_domain_accessors();
    test_domain_sind_vind();
    test_domain_get_cell_index();
    test_domain_to_cell_coordinates();
    test_domain_pos_vind_pos_sind();
    test_domain_is_inside_node_linear();
    test_domain_is_inside_node_periodic();

    test::print_summary();

    return test::get_failed();
}
