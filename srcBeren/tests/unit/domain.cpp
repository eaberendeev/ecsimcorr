#include "grid/Domain.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "core/Vector3.h"
#include "test_case.h"

#if defined(NLOHMANN_JSON_VERSION_MAJOR) || defined(EIGEN_WORLD_VERSION)
#error "Grid/domain headers must not import JSON or Eigen"
#endif

// Grid indices include ghost layers. Yee E components have one half-cell
// shift; B components have two. Geometric examples below use known points
// and mirror images, independently of particle storage and field solvers.
// =============================================================================
// Тесты для Vector3R
// =============================================================================
void test_vector3_construction() {
    Vector3R v1;
    check::equal(v1.x(), 0.0, 1e-10, "Vector3R default constructor x");
    check::equal(v1.y(), 0.0, 1e-10, "Vector3R default constructor y");
    check::equal(v1.z(), 0.0, 1e-10, "Vector3R default constructor z");

    Vector3R v2(1.0, 2.0, 3.0);
    check::equal(v2.x(), 1.0, 1e-10, "Vector3R param constructor x");
    check::equal(v2.y(), 2.0, 1e-10, "Vector3R param constructor y");
    check::equal(v2.z(), 3.0, 1e-10, "Vector3R param constructor z");

    Vector3R v3(5.0);
    check::equal(v3.x(), 5.0, 1e-10, "Vector3R uniform constructor x");
    check::equal(v3.y(), 5.0, 1e-10, "Vector3R uniform constructor y");
    check::equal(v3.z(), 5.0, 1e-10, "Vector3R uniform constructor z");
}

void test_vector3_operations() {
    Vector3R v1(1.0, 2.0, 3.0);
    Vector3R v2(4.0, 5.0, 6.0);

    Vector3R sum = v1 + v2;
    check::equal(sum.x(), 5.0, 1e-10, "Vector3 addition x");
    check::equal(sum.y(), 7.0, 1e-10, "Vector3 addition y");
    check::equal(sum.z(), 9.0, 1e-10, "Vector3 addition z");

    Vector3R diff = v2 - v1;
    check::equal(diff.x(), 3.0, 1e-10, "Vector3 subtraction x");
    check::equal(diff.y(), 3.0, 1e-10, "Vector3 subtraction y");
    check::equal(diff.z(), 3.0, 1e-10, "Vector3 subtraction z");

    Vector3R scaled = v1 * 2.0;
    check::equal(scaled.x(), 2.0, 1e-10, "Vector3 scalar mult x");
    check::equal(scaled.y(), 4.0, 1e-10, "Vector3 scalar mult y");
    check::equal(scaled.z(), 6.0, 1e-10, "Vector3 scalar mult z");

    double dot = v1.dot(v2);
    check::equal(dot, 32.0, 1e-10, "Vector3 dot product");

    Vector3R cross = v1.cross(v2);
    check::equal(cross.x(), -3.0, 1e-10, "Vector3 cross product x");
    check::equal(cross.y(), 6.0, 1e-10, "Vector3 cross product y");
    check::equal(cross.z(), -3.0, 1e-10, "Vector3 cross product z");
}

// =============================================================================
// Тесты для Grid
// =============================================================================
void test_grid_construction() {
    Vector3R cell_size(0.5, 0.5, 0.5);
    Vector3I num_cells(10, 10, 10);
    int ghost_cells = 2;

    Grid grid(cell_size, num_cells, ghost_cells);

    check::equal(grid.cell_size().x(), 0.5, 1e-10, "Grid cell_size x");
    check::equal(grid.cell_size().y(), 0.5, 1e-10, "Grid cell_size y");
    check::equal(grid.cell_size().z(), 0.5, 1e-10, "Grid cell_size z");

    check::equal(grid.num_cells().x(), 10, "Grid num_cells x");
    check::equal(grid.num_cells().y(), 10, "Grid num_cells y");
    check::equal(grid.num_cells().z(), 10, "Grid num_cells z");

    check::equal(grid.ghost_cells(), 2, "Grid ghost_cells");

    // size = num_cells + 2*ghost + 1 = 10 + 4 + 1 = 15
    check::equal(grid.size().x(), 15, "Grid size x");
    check::equal(grid.size().y(), 15, "Grid size y");
    check::equal(grid.size().z(), 15, "Grid size z");

    double expected_volume = 0.5 * 0.5 * 0.5;
    check::equal(grid.cell_volume(), expected_volume, 1e-10, "Grid cell_volume");
}

void test_grid_indexing() {
    Vector3R cell_size(1.0, 1.0, 1.0);
    Vector3I num_cells(3, 4, 5);
    int ghost_cells = 1;
    // size 6 7 8
    Grid grid(cell_size, num_cells, ghost_cells);

    // sind(i, j, k) = i * size.y * size.z + j * size.z + k
    // 2, 3, 1
    int expected_sind = 2 * 7 * 8 + 3 * 8 + 1;
    check::equal(grid.sind(2, 3, 1), expected_sind, "Grid sind");

    // vind(i, j, k, d, nd) = d + nd * (i * size.y * size.z + j * size.z + k)
    int expected_vind = 2 + 3 * (1 * 7 * 8 + 2 * 8 + 3);
    check::equal(grid.vind(1, 2, 3, 2), expected_vind, "Grid vind");
}

void test_grid_coordinates() {
    Vector3R cell_size(0.5, 0.5, 0.5);
    Vector3I num_cells(11, 12, 13);
    int ghost_cells = 1;

    Grid grid(cell_size, num_cells, ghost_cells);

    Vector3R world(1.0, 2.0, 3.0);
    Vector3R cell = grid.to_cell_coordinates(world);
    check::equal(cell.x(), 2.0, 1e-10, "Grid to_cell_coordinates x");
    check::equal(cell.y(), 4.0, 1e-10, "Grid to_cell_coordinates y");
    check::equal(cell.z(), 6.0, 1e-10, "Grid to_cell_coordinates z");

    Vector3R coord(1.25, 1.75, 2.0);
    Vector3I idx = grid.get_cell_index(coord);
    // idx = int(coord / cell_size + ghost) = int(2.5 + 1) = 3
    // idy = int(coord / cell_size + ghost) = int(3.5 + 1) = 4
    // idz = int(coord / cell_size + ghost) = int(4 + 1) = 5
    check::equal(idx.x(), 3, "Grid get_cell_index x");
    check::equal(idx.y(), 4, "Grid get_cell_index y");
    check::equal(idx.z(), 5, "Grid get_cell_index z");
}

// =============================================================================
// Тесты для Geometry
// =============================================================================
void test_geometry_box() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(11.0, 13.0, 12.0);

    Geometry geom(box_min, box_max, false);

    Vector3R inside(5.0, 5.0, 5.0);
    check::require(geom.contains(inside, 1e-10), "Geometry contains point inside box");

    Vector3R outside(-0.01, 3.0, 12.0);
    check::require(!geom.contains(outside, 1e-10), "Geometry rejects point outside box (x < min)");

    Vector3R outside2(11.0, 5.0, 5.0);
    check::require(!geom.contains(outside2), "Geometry rejects point on boundary (x = max)");

    Vector3R on_min(0.001, 0.0001, 0.00001);
    check::require(geom.contains(on_min, 1e-10), "Geometry contains point near min boundary");
}

void test_geometry_cylinder() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(13.0, 13.0, 12.0);

    Geometry geom(box_min, box_max, true);
    geom.cyl_center = Vector3R(6.5, 6.5, 0.0);
    geom.cyl_radius = 6.5;

    // Точка внутри цилиндра (на оси)
    Vector3R inside_cyl(5.0, 5.0, 5.0);
    check::require(geom.contains(inside_cyl, 1e-10), "Geometry contains point inside cylinder");

    // Точка вне цилиндра (слишком далеко от оси)
    Vector3R outside_cyl(-1.5, 1.5, 5.0);
    check::require(!geom.contains(outside_cyl, 1e-10), "Geometry rejects point outside cylinder");

    // Точка на границе цилиндра
    Vector3R on_boundary(0.0, 0.0, 1.0);
    check::require(!geom.contains(on_boundary, 1e-10), "Geometry rejects point on cylinder boundary");

    // Точка внутри box, но вне цилиндра по x
    Vector3R outside_box(1.5, 1.5, 5.0);
    // distance^2 = (1.5-6.5)^2 + (1.5-6.5)^2 = 25 + 25 = 50 > 6.5*6.5
    check::require(!geom.contains(outside_box, 1e-10), "Geometry rejects point outside box");
}

void test_geometry_faces() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(8.0, 13.0, 20.0);

    Geometry geom(box_min, box_max, false);

    Vector3R p(8.0, 3.0, 6.0);
    check::require(!geom.is_outside_face(Face::XMIN, p, 1e-10), "Point not outside XMIN");

    Vector3R p_xmin(-0.01, 3.0, 15.0);
    check::require(geom.is_outside_face(Face::XMIN, p_xmin, 1e-10), "Point outside XMIN");

    Vector3R p_xmax(8.0, 2.0, 0.01);
    check::require(geom.is_outside_face(Face::XMAX, p_xmax, 1e-10), "Point outside XMAX");

    Vector3R p_ymin(5.0, -0.00004, 5.0);
    check::require(geom.is_outside_face(Face::YMIN, p_ymin, 1e-10), "Point outside YMIN");

    Vector3R p_ymax(-1.0, 13.0, 5.0);
    check::require(geom.is_outside_face(Face::YMAX, p_ymax, 1e-10), "Point outside YMAX");

    Vector3R p_zmin(-5.0, 5.0, -0.0001);
    check::require(geom.is_outside_face(Face::ZMIN, p_zmin, 1e-10), "Point outside ZMIN");

    Vector3R p_zmax(5.0, 5.0, 20.0000000001);
    check::require(geom.is_outside_face(Face::ZMAX, p_zmax, 1e-10), "Point outside ZMAX");
}

// =============================================================================
// Тесты для Domain
// =============================================================================
void test_domain_init() {
    Domain domain;
    Vector3I num_cells(105, 10, 5);
    Vector3R cell_size(0.1, 0.2, 0.3);
    int ghost = 3;

    domain.init(num_cells, cell_size, ghost);

    check::equal(domain.grid.cell_size().x(), 0.1, 1e-10, "Domain grid cell_size x");
    check::equal(domain.grid.num_cells().x(), 105, "Domain grid num_cells x");
    check::equal(domain.grid.ghost_cells(), 3, "Domain grid ghost_cells");

    // box_max = num_cells * cell_size
    check::equal(domain.geom.box_max.x(), 10.5, 1e-10, "Domain geom box_max x");
    check::equal(domain.geom.box_max.y(), 2.0, 1e-10, "Domain geom box_max y");
    check::equal(domain.geom.box_max.z(), 1.5, 1e-10, "Domain geom box_max z");
}

void test_domain_cylinder() {
    Domain domain;
    Vector3I num_cells(20, 20, 20);
    Vector3R cell_size(0.2, 0.2, 0.2);

    domain.init(num_cells, cell_size, 1);
    domain.set_cylinder(Vector3R(2.0, 2.0, 2.0), 2.0);

    check::require(domain.geom.use_cylinder, "Domain cylinder enabled");
    check::equal(domain.geom.cyl_center.x(), 2.0, 1e-10, "Domain cyl_center x");
    check::equal(domain.geom.cyl_center.y(), 2.0, 1e-10, "Domain cyl_center y");
    check::equal(domain.geom.cyl_radius, 2.0, 1e-10, "Domain cyl_radius");

    // Точка внутри цилиндра
    Vector3R inside(2.5, 1.5, 1.5);
    check::require(domain.contains(inside, 1e-10), "Domain contains point inside cylinder");

    // Точка вне цилиндра
    Vector3R outside(3.8, 3.5, 3.0);
    check::require(!domain.contains(outside, 1e-10), "Domain rejects point outside cylinder");
}

void test_domain_node_position() {
    Domain domain;
    Vector3I num_cells(12, 13, 11);
    Vector3R cell_size(0.5, 0.5, 0.5);

    domain.init(num_cells, cell_size, 1);

    // Тест для электрического поля (Ex компонент)
    // Ex сдвинута на +dx/2 по X
    Vector3R pos_ex = domain.get_node_position(1, 1, 1, FieldType::ELECTRIC, 0);
    // dx=0.5: x=(i-ghost+0.5)*dx=0.25; y=z=0.
    check::equal(pos_ex.x(), 0.25, 1e-10, "Domain Ex node position x");
    check::equal(pos_ex.y(), 0.0, 1e-10, "Domain Ex node position y");
    check::equal(pos_ex.z(), 0.0, 1e-10, "Domain Ex node position z");

    Vector3R pos_ey = domain.get_node_position(0, 2, 1, FieldType::ELECTRIC, 1);
    check::equal(pos_ey.x(), -0.5, 1e-10, "Domain Ex node position x");
    check::equal(pos_ey.y(), 0.75, 1e-10, "Domain Ex node position y");
    check::equal(pos_ey.z(), 0.0, 1e-10, "Domain Ex node position z");

    Vector3R pos_ez = domain.get_node_position(0, 1, 1, FieldType::ELECTRIC, 2);
    check::equal(pos_ez.x(), -0.5, 1e-10, "Domain Ex node position x");
    check::equal(pos_ez.y(), 0.0, 1e-10, "Domain Ex node position y");
    check::equal(pos_ez.z(), 0.25, 1e-10, "Domain Ex node position z");

    // Тест для магнитного поля (Bx компонент)
    // Bx сдвинута на +dy/2 и +dz/2
    Vector3R pos_bx = domain.get_node_position(1, 1, 1, FieldType::MAGNETIC, 0);
    check::equal(pos_bx.x(), 0.0, 1e-10, "Domain Bx node position x");
    check::equal(pos_bx.y(), 0.25, 1e-10, "Domain Bx node position y");
    check::equal(pos_bx.z(), 0.25, 1e-10, "Domain Bx node position z");
    Vector3R pos_by = domain.get_node_position(1, 1, 1, FieldType::MAGNETIC, 1);
    check::equal(pos_by.x(), 0.25, 1e-10, "Domain Bx node position x");
    check::equal(pos_by.y(), 0.0, 1e-10, "Domain Bx node position y");
    check::equal(pos_by.z(), 0.25, 1e-10, "Domain Bx node position z");
    Vector3R pos_bz = domain.get_node_position(1, 1, 1, FieldType::MAGNETIC, 2);
    check::equal(pos_bz.x(), 0.25, 1e-10, "Domain Bx node position x");
    check::equal(pos_bz.y(), 0.25, 1e-10, "Domain Bx node position y");
    check::equal(pos_bz.z(), 0.0, 1e-10, "Domain Bx node position z");
    // Плотность задана в узлах без Yee-сдвига.
    Vector3R pos_density = domain.get_node_position(1, 1, 1, FieldType::DENSITY, 0);
    check::equal(pos_density.x(), 0.0, 1e-10, "Domain density node position x");
    check::equal(pos_density.y(), 0.0, 1e-10, "Domain density node position y");
    check::equal(pos_density.z(), 0.0, 1e-10, "Domain density node position z");
}

void test_domain_is_inside_node() {
    Domain domain;
    Vector3I num_cells(14, 14, 7);
    Vector3R cell_size(0.5, 0.5, 0.5);

    domain.init(num_cells, cell_size, 1);
    auto size = domain.grid.size();

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                bool inside_ex = domain.is_inside_node(i, j, k, FieldType::ELECTRIC, 0);
                if (i < 1 || j < 2 || k < 2 || i > size.x() - 3 || j > size.y() - 3 || k > size.z() - 3) {
                    check::require(!inside_ex, "Domain node outside box");
                } else {
                    check::require(inside_ex, "Domain node inside box");
                }
                bool inside_ey = domain.is_inside_node(i, j, k, FieldType::ELECTRIC, 1);
                if (i < 2 || j < 1 || k < 2 || i > size.x() - 3 || j > size.y() - 3 || k > size.z() - 3) {
                    check::require(!inside_ey, "Domain node outside box");
                } else {
                    check::require(inside_ey, "Domain node inside box");
                }
                bool inside_ez = domain.is_inside_node(i, j, k, FieldType::ELECTRIC, 2);
                if (i < 2 || j < 2 || k < 1 || i > size.x() - 3 || j > size.y() - 3 || k > size.z() - 3) {
                    check::require(!inside_ez, "Domain node outside box");
                } else {
                    check::require(inside_ez, "Domain node inside box");
                }
                bool inside_bx = domain.is_inside_node(i, j, k, FieldType::MAGNETIC, 0);
                if (i < 2 || j < 1 || k < 1 || i > size.x() - 3 || j > size.y() - 3 || k > size.z() - 3) {
                    check::require(!inside_bx, "Domain node outside box");
                } else {
                    check::require(inside_bx, "Domain node inside box");
                }
                bool inside_by = domain.is_inside_node(i, j, k, FieldType::MAGNETIC, 1);
                if (i < 1 || j < 2 || k < 1 || i > size.x() - 3 || j > size.y() - 3 || k > size.z() - 3) {
                    check::require(!inside_by, "Domain node outside box");
                } else {
                    check::require(inside_by, "Domain node inside box");
                }
                bool inside_bz = domain.is_inside_node(i, j, k, FieldType::MAGNETIC, 2);
                if (i < 1 || j < 1 || k < 2 || i > size.x() - 3 || j > size.y() - 3 || k > size.z() - 3) {
                    check::require(!inside_bz, "Domain node outside box");
                } else {
                    check::require(inside_bz, "Domain node inside box");
                }
            }
        }
    }

    domain.set_cylinder(Vector3R(3.5, 3.5, 0.0), 3.5);

    // Точка внутри цилиндра
    bool inside = domain.is_inside_node(5, 6, 2, FieldType::DENSITY, 0);
    check::require(inside, "Domain node inside cylinder");

    // Точка вне цилиндра
    bool outside = domain.is_inside_node(2, 3, 3, FieldType::DENSITY, 0);
    check::require(!outside, "Domain node outside cylinder");
}

// =============================================================================
// Дополнительные тесты для Grid
// =============================================================================
void test_grid_origin_and_total() {
    Vector3R cell_size(1.0, 2.0, 3.0);
    Vector3I num_cells(5, 6, 7);
    int ghost_cells = 2;

    Grid grid(cell_size, num_cells, ghost_cells);

    check::equal(grid.origin().x(), 0.0, 1e-10, "Grid origin x");
    check::equal(grid.origin().y(), 0.0, 1e-10, "Grid origin y");
    check::equal(grid.origin().z(), 0.0, 1e-10, "Grid origin z");

    // size = (5+4+1, 6+4+1, 7+4+1) = (10, 11, 12)
    int expected_total = 10 * 11 * 12;
    check::equal(grid.total_size(), expected_total, "Grid total_size");
}

void test_grid_pos_vind() {
    Vector3R cell_size(1.0, 1.0, 1.0);
    Vector3I num_cells(3, 4, 5);
    int ghost_cells = 1;
    // size: 6, 7, 8; dims_ = {6, 7, 8, 3}

    Grid grid(cell_size, num_cells, ghost_cells);

    int index = 0;
    check::equal(grid.pos_vind(index, 0), 0, "Grid pos_vind index=0 n=0");
    check::equal(grid.pos_vind(index, 1), 0, "Grid pos_vind index=0 n=1");
    check::equal(grid.pos_vind(index, 2), 0, "Grid pos_vind index=0 n=2");
    check::equal(grid.pos_vind(index, 3), 0, "Grid pos_vind index=0 n=3");

    // index = d + nd*(i*sy*sz + j*sz + k) = 2 + 3*(1*7*8 + 2*8 + 3) = 2 + 3*(56+16+3) = 2 + 225 = 227
    index = 227;
    check::equal(grid.pos_vind(index, 0), 1, "Grid pos_vind index=227 n=0 (i)");
    check::equal(grid.pos_vind(index, 1), 2, "Grid pos_vind index=227 n=1 (j)");
    check::equal(grid.pos_vind(index, 2), 3, "Grid pos_vind index=227 n=2 (k)");
    check::equal(grid.pos_vind(index, 3), 2, "Grid pos_vind index=227 n=3 (d)");
}

void test_grid_pos_sind() {
    Vector3R cell_size(1.0, 1.0, 1.0);
    Vector3I num_cells(3, 4, 5);
    int ghost_cells = 1;
    // size: 6, 7, 8

    Grid grid(cell_size, num_cells, ghost_cells);

    int index = 0;
    check::equal(grid.pos_sind(index, 0), 0, "Grid pos_sind index=0 n=0");
    check::equal(grid.pos_sind(index, 1), 0, "Grid pos_sind index=0 n=1");
    check::equal(grid.pos_sind(index, 2), 0, "Grid pos_sind index=0 n=2");

    // index = i*sy*sz + j*sz + k = 1*7*8 + 2*8 + 3 = 56 + 16 + 3 = 75
    index = 75;
    check::equal(grid.pos_sind(index, 0), 1, "Grid pos_sind index=75 n=0 (i)");
    check::equal(grid.pos_sind(index, 1), 2, "Grid pos_sind index=75 n=1 (j)");
    check::equal(grid.pos_sind(index, 2), 3, "Grid pos_sind index=75 n=2 (k)");
}

void test_grid_get_field_node_index() {
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
    check::equal(idx_ex.x(), 2, "Grid Ex node index x");
    check::equal(idx_ex.y(), 5, "Grid Ex node index y");
    check::equal(idx_ex.z(), 7, "Grid Ex node index z");

    // Bx: shift = (0, 0.5, 0.5) in cell units
    // x = 2.0 - 0 = 2.0, floor = 2, i = 3
    // y = 4.0 - 0.5 = 3.5, floor = 3, j = 4
    // z = 6.0 - 0.5 = 5.5, floor = 5, k = 6
    Vector3I idx_bx = grid.get_field_node_index(world, FieldType::MAGNETIC, 0);
    check::equal(idx_bx.x(), 3, "Grid Bx node index x");
    check::equal(idx_bx.y(), 4, "Grid Bx node index y");
    check::equal(idx_bx.z(), 6, "Grid Bx node index z");

    // Density: shift = (0, 0, 0)
    // x = 2.0, floor = 2, i = 3
    // y = 4.0, floor = 4, j = 5
    // z = 6.0, floor = 6, k = 7
    Vector3I idx_rho = grid.get_field_node_index(world, FieldType::DENSITY, 0);
    check::equal(idx_rho.x(), 3, "Grid density node index x");
    check::equal(idx_rho.y(), 5, "Grid density node index y");
    check::equal(idx_rho.z(), 7, "Grid density node index z");
}

// =============================================================================
// Дополнительные тесты для Geometry
// =============================================================================
void test_geometry_in_cylinder() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);

    Geometry geom;
    geom.init(box_min, box_max, false);
    check::require(!geom.in_cylinder(Vector3R(5.0, 5.0, 5.0), 1e-10), "in_cylinder false when use_cylinder=false");

    geom.use_cylinder = true;
    geom.cyl_center = Vector3R(5.0, 5.0, 0.0);
    geom.cyl_radius = 3.0;

    Vector3R inside(5.0, 5.0, 5.0);
    check::require(geom.in_cylinder(inside, 1e-10), "in_cylinder at center");

    Vector3R near_boundary(5.0, 7.9, 5.0);
    check::require(geom.in_cylinder(near_boundary, 1e-10), "in_cylinder near boundary");

    Vector3R outside(5.0, 8.1, 5.0);
    check::require(!geom.in_cylinder(outside, 1e-10), "in_cylinder outside");
}

void test_geometry_contains_periodic() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    bool periodic[3] = {false, false, false};
    // x outside
    check::require(!geom.contains(Vector3R(-1.0, 5.0, 5.0), periodic, 1e-10), "contains periodic: x < min rejected");
    // periodic x — пропускает границу
    bool periodic_x[3] = {true, false, false};
    check::require(geom.contains(Vector3R(-1.0, 5.0, 5.0), periodic_x, 1e-10),
                   "contains periodic: x periodic skips x min");
    // periodic all — любая точка внутри
    bool periodic_all[3] = {true, true, true};
    check::require(geom.contains(Vector3R(-999.0, 999.0, -999.0), periodic_all, 1e-10),
                   "contains periodic: all periodic accepts anything");
}

void test_geometry_contains_ignoring_face() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Точка вне XMIN, но XMIN игнорируется — внутри
    Vector3R p(-1.0, 5.0, 5.0);
    check::require(geom.contains_ignoring_face(Face::XMIN, p, 1e-10), "contains_ignoring_face XMIN");
    check::require(!geom.contains_ignoring_face(Face::XMAX, p, 1e-10), "contains_ignoring_face XMAX rejects x < min");

    // Точка вне XMAX, но XMAX игнорируется — внутри
    Vector3R p2(10.5, 5.0, 5.0);
    check::require(geom.contains_ignoring_face(Face::XMAX, p2, 1e-10), "contains_ignoring_face XMAX");

    // CYLINDER
    geom.use_cylinder = true;
    geom.cyl_center = Vector3R(5.0, 5.0, 0.0);
    geom.cyl_radius = 3.0;
    Vector3R outside_cyl(5.0, 9.0, 5.0);
    check::require(!geom.contains_ignoring_face(Face::XMIN, outside_cyl, 1e-10),
                   "contains_ignoring_face rejects outside cylinder");
    check::require(geom.contains_ignoring_face(Face::CYLINDER, outside_cyl, 1e-10),
                   "contains_ignoring_face CYLINDER ignores cylinder check");
}

void test_geometry_is_outside_only_face() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Точка чуть левее XMIN, все остальные грани внутри
    Vector3R p(-0.1, 5.0, 5.0);
    check::require(geom.is_outside_only_face(Face::XMIN, p, 1e-10), "is_outside_only_face XMIN only");

    // Точка вне XMIN и YMIN одновременно
    Vector3R p2(-0.1, -0.1, 5.0);
    check::require(!geom.is_outside_only_face(Face::XMIN, p2, 1e-10),
                   "is_outside_only_face: outside two faces, not only XMIN");
}

void test_geometry_reflect_from_face() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Отражаем от XMIN: x' = 2*0 - (-1) = 1
    Vector3R p(-1.0, 5.0, 5.0);
    Vector3R refl = geom.reflect_from_face(Face::XMIN, p);
    check::equal(refl.x(), 1.0, 1e-10, "reflect_from_face XMIN x");
    check::equal(refl.y(), 5.0, 1e-10, "reflect_from_face XMIN y");
    check::equal(refl.z(), 5.0, 1e-10, "reflect_from_face XMIN z");

    // Отражаем от XMAX: x' = 2*10 - 11 = 9
    Vector3R p2(11.0, 3.0, 7.0);
    Vector3R refl2 = geom.reflect_from_face(Face::XMAX, p2);
    check::equal(refl2.x(), 9.0, 1e-10, "reflect_from_face XMAX x");

    // YMIN
    Vector3R p3(5.0, -2.0, 5.0);
    Vector3R refl3 = geom.reflect_from_face(Face::YMIN, p3);
    check::equal(refl3.y(), 2.0, 1e-10, "reflect_from_face YMIN y");

    // YMAX
    Vector3R p4(5.0, 12.0, 5.0);
    Vector3R refl4 = geom.reflect_from_face(Face::YMAX, p4);
    check::equal(refl4.y(), 8.0, 1e-10, "reflect_from_face YMAX y");

    // ZMIN
    Vector3R p5(5.0, 5.0, -3.0);
    Vector3R refl5 = geom.reflect_from_face(Face::ZMIN, p5);
    check::equal(refl5.z(), 3.0, 1e-10, "reflect_from_face ZMIN z");

    // ZMAX
    Vector3R p6(5.0, 5.0, 13.0);
    Vector3R refl6 = geom.reflect_from_face(Face::ZMAX, p6);
    check::equal(refl6.z(), 7.0, 1e-10, "reflect_from_face ZMAX z");
}

void test_geometry_reflect_cylinder() {
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
    check::equal(refl.x(), 8.0, 1e-10, "reflect CYLINDER x");
    check::equal(refl.y(), 5.0, 1e-10, "reflect CYLINDER y");
    check::equal(refl.z(), 3.0, 1e-10, "reflect CYLINDER z (unchanged)");
}

void test_geometry_reflect_from_boundary() {
    Vector3R box_min(0.0, 0.0, 0.0);
    Vector3R box_max(10.0, 10.0, 10.0);
    Geometry geom(box_min, box_max, false);

    // Точка за XMIN
    Vector3R p(-1.0, 5.0, 5.0);
    Vector3R refl = geom.reflect_from_boundary(p);
    check::equal(refl.x(), 1.0, 1e-10, "reflect_boundary xmin x");
    check::require(geom.contains(refl, 1e-12), "reflect_boundary: result inside box");

    // Точка за двумя гранями (XMIN и YMAX)
    Vector3R p2(-2.0, 13.0, 5.0);
    Vector3R refl2 = geom.reflect_from_boundary(p2);
    check::require(geom.contains(refl2, 1e-12), "reflect_boundary: result inside after double reflect");

    // Точка за XMAX, ZMIN
    Vector3R p3(12.0, 5.0, -1.0);
    Vector3R refl3 = geom.reflect_from_boundary(p3);
    check::require(geom.contains(refl3, 1e-12), "reflect_boundary: triple axis check inside");
}

// =============================================================================
// Дополнительные тесты для Domain
// =============================================================================
void test_domain_is_inside_node_linear() {
    Domain domain;
    Vector3I num_cells(5, 5, 5);
    Vector3R cell_size(1.0, 1.0, 1.0);

    domain.init(num_cells, cell_size, 1);

    // vind(3, 3, 3, 0) для ELECTRIC (Ex) — внутри бокса
    int idx_ex = domain.vind(3, 3, 3, 0);
    check::require(domain.is_inside_node(idx_ex, FieldType::ELECTRIC), "Domain is_inside_node linear Ex inside");
    // vind(0, 3, 3, 0) — снаружи
    int idx_out = domain.vind(0, 3, 3, 0);
    check::require(!domain.is_inside_node(idx_out, FieldType::ELECTRIC), "Domain is_inside_node linear Ex outside");
    // vind(3, 3, 3, 1) для ELECTRIC (Ey) — внутри
    int idx_ey = domain.vind(3, 3, 3, 1);
    check::require(domain.is_inside_node(idx_ey, FieldType::ELECTRIC), "Domain is_inside_node linear Ey inside");
}

void test_domain_is_inside_node_periodic() {
    Domain domain;
    Vector3I num_cells(5, 5, 5);
    Vector3R cell_size(1.0, 1.0, 1.0);

    domain.init(num_cells, cell_size, 1);

    bool no_periodic[3] = {false, false, false};
    bool periodic_x[3] = {true, false, false};

    // Крайний узел Ex (i=0) — за левой границей без periodicity
    check::require(!domain.is_inside_node_periodic(0, 3, 3, FieldType::ELECTRIC, 0, no_periodic),
                   "is_inside_node_periodic: Ex at i=0 non-periodic rejected");
    // С periodicity по X — принимается
    check::require(domain.is_inside_node_periodic(0, 3, 3, FieldType::ELECTRIC, 0, periodic_x),
                   "is_inside_node_periodic: Ex at i=0 periodic x accepted");
}

// =============================================================================
// Главная функция
// =============================================================================
void noncubic_layout_roundtrips_and_scalar_nodes() {
    const Vector3I cells(3, 4, 5);
    for (int ghost : {0, 1, 2}) {
        Domain domain;
        domain.init(cells, Vector3R(0.5, 1, 2), ghost);
        const Vector3I sizes = domain.size();
        for (int i = 0; i < sizes.x(); ++i)
            for (int j = 0; j < sizes.y(); ++j)
                for (int k = 0; k < sizes.z(); ++k) {
                    const int scalar = i * sizes.y() * sizes.z() + j * sizes.z() + k;
                    check::equal(domain.sind(i, j, k), scalar, "scalar row-major index");
                    check::equal(domain.pos_sind(scalar, 0), i, "scalar x coordinate");
                    check::equal(domain.pos_sind(scalar, 1), j, "scalar y coordinate");
                    check::equal(domain.pos_sind(scalar, 2), k, "scalar z coordinate");
                    // Scalar density nodes at box faces are outside the active
                    // interior selected by the existing positive epsilon.
                    const bool interior = i > ghost && i < ghost + cells.x() && j > ghost && j < ghost + cells.y() &&
                                          k > ghost && k < ghost + cells.z();
                    check::require(domain.is_inside_node(scalar, FieldType::DENSITY) == interior,
                                   "density index was decoded using vector components");
                    for (int d = 0; d < 3; ++d) {
                        const int vector = 3 * scalar + d;
                        check::equal(domain.vind(i, j, k, d), vector, "vector row-major index");
                        const std::array<int, 4> decoded = domain.grid.pos_vind_range(vector);
                        check::require(decoded == std::array<int, 4>{i, j, k, d}, "vector index round trip");
                        for (FieldType field :
                             {FieldType::ELECTRIC, FieldType::MAGNETIC, FieldType::DENSITY, FieldType::CURRENT}) {
                            const Vector3R position = domain.get_node_position(i, j, k, field, d);
                            const Vector3I node = domain.grid.get_field_node_index(position, field, d);
                            check::require(node.x() == i && node.y() == j && node.z() == k,
                                           "Yee coordinate/index round trip changed with ghost width");
                        }
                    }
                }
    }
}

int main() {
    return check::run([] {
        check::scenario("Noncubic scalar/vector layouts, Yee coordinates and 0/1/2 ghost layers",
                        noncubic_layout_roundtrips_and_scalar_nodes);
        check::scenario("Vector construction and vector algebra", [] {
            test_vector3_construction();
            test_vector3_operations();
        });
        check::scenario("Grid sizes, cell coordinates and ghost offsets", [] {
            test_grid_construction();
            test_grid_indexing();
            test_grid_coordinates();
            test_grid_origin_and_total();
        });
        check::scenario("Flattened indices and Yee node coordinates", [] {
            test_grid_pos_vind();
            test_grid_pos_sind();
            test_grid_get_field_node_index();
            test_domain_node_position();
        });
        check::scenario("Box and cylinder containment, including periodic faces", [] {
            test_geometry_box();
            test_geometry_cylinder();
            test_geometry_faces();
            test_geometry_in_cylinder();
            test_geometry_contains_periodic();
            test_geometry_contains_ignoring_face();
            test_geometry_is_outside_only_face();
        });
        check::scenario("Geometric reflection from planes and cylinder", [] {
            test_geometry_reflect_from_face();
            test_geometry_reflect_cylinder();
            test_geometry_reflect_from_boundary();
        });
        check::scenario("Domain geometry and active Yee nodes", [] {
            test_domain_init();
            test_domain_cylinder();
            test_domain_is_inside_node();
            test_domain_is_inside_node_linear();
            test_domain_is_inside_node_periodic();
        });
    });
}
