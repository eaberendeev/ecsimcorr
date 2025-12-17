#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "containers.h"
#include "service.h"
#include "voxel_traversal.h"
#include "interpolation.h"
void ParticlesArray::move(double dt) {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.move(dt);
        }
    }
}

void ParticlesArray::predict_velocity(const Field3d& fieldE,
                                      const Field3d& fieldEp,
                                      const Field3d& fieldB,
                                      const Domain& domain, const double dt,
                                      ShapeType type) {
    if(is_neutral()) return;

    switch (type) {
        case ShapeType::NGP:
            predict_velocity_impl_ngp(fieldE, fieldEp, fieldB, domain, dt);
            break;
        case ShapeType::Linear:
            predict_velocity_impl_linear(fieldE, fieldEp, fieldB, domain, dt);
            break;
        case ShapeType::Quadratic:
            std::cout
                << "Predict velocity for quadratic shape not implemented\n";
            exit(-1);
    }
}

void ParticlesArray::predict_velocity_impl_linear(const Field3d& fieldE,
                                      const Field3d& fieldEp,
                                      const Field3d& fieldB,
                                      const Domain& domain, const double dt) {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            const auto coord = particle.coord;
            const auto velocity = particle.velocity;
            Vector3R E = interpolateE_linear(fieldE, to_cell_coordinates(coord));
            const Vector3R B =
                interpolateB_linear(fieldB, to_cell_coordinates(coord));
            const Vector3R En =
                interpolateE_linear(fieldEp, to_cell_coordinates(coord));
            E = 0.5 * (E + En);
            const auto beta = dt * charge / _mass;
            const auto alpha = 0.5 * beta * B.norm();
            const auto alpha2 = alpha * alpha;
            const auto h = B.normalized();

            const auto v12 =
                (1. / (1. + alpha2)) *
                (velocity + alpha * velocity.cross(h) + alpha2 * h.dot(velocity) * h +
                     0.5 * beta *
                         (E + alpha * E.cross(h) + alpha2 * E.dot(h) * h));

            particle.velocity = 2. * v12 - velocity;
        }
    }
}

void ParticlesArray::predict_velocity_impl_ngp(
    const Field3d& fieldE, const Field3d& fieldEp, const Field3d& fieldB,
    [[maybe_unused]] const Domain& domain, const double dt) {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            const auto coord = particle.coord;
            const auto velocity = particle.velocity;
            const auto normalized_coord = to_cell_coordinates(coord);
            Vector3R E = interpolateE_ngp(fieldE, normalized_coord);
            const Vector3R B = interpolateB_ngp(fieldB, normalized_coord);
            const Vector3R En = interpolateE_ngp(fieldEp, normalized_coord);
            E = 0.5 * (E + En);
            const auto beta = dt * charge / _mass;
            const auto alpha = 0.5 * beta * B.norm();
            const auto alpha2 = alpha * alpha;
            const auto h = B.normalized();

            const auto v12 =
                (1. / (1. + alpha2)) *
                (velocity + alpha * velocity.cross(h) + alpha2 * h.dot(velocity) * h +
                     0.5 * beta *
                         (E + alpha * E.cross(h) + alpha2 * E.dot(h) * h));

            particle.velocity = 2. * v12 - velocity;
        }
    }
}

// Template specializations need to be explicitly instantiated in the cpp file
template void ParticlesArray::move_and_calc_current_impl<Shape, 2>(
    const double dt, Field3d& fieldJ);
template void ParticlesArray::move_and_calc_current_impl<Shape2, 2>(
    const double dt, Field3d& fieldJ);

void ParticlesArray::move_and_calc_current(const double dt, Field3d& fieldJ,
                                           ShapeType type) {
    if (is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            std::cout << "Move and calc current for NGP is not supported\n" << std::endl;
            exit(-1);
        case ShapeType::Linear:
            move_and_calc_current_impl<Shape, 2>(dt, fieldJ);
            break;
        case ShapeType::Quadratic:
            move_and_calc_current_impl<Shape2, 2>(dt, fieldJ);
            break;
    }
}

template <ParticlesArray::ShapeFunction ShapeFn, int ShapeSize>
void ParticlesArray::move_and_calc_current_impl(const double dt,
                                                Field3d& fieldJ) {
    constexpr auto SMAX = 2 * ShapeSize;

    const double qx = charge * xCellSize / (6 * dt) * _mpw;
    const double qy = charge * yCellSize / (6 * dt) * _mpw;
    const double qz = charge * zCellSize / (6 * dt) * _mpw;

// TODO: change base_ to cell index from ParticlesData
#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < size(); ++pk) {
        auto& particles = particlesData(pk);
        if (particles.empty()) {
            continue;
        }

        ParticleShape<ShapeFn, SMAX> start_shape;
        ParticleShape<ShapeFn, SMAX> end_shape;
        CurrentBuffer<SMAX> curBuf, cellBuf;
        cellBuf.zero();
        start_shape.fill_zero();

        for (auto& particle : particles) {
            Vector3R start = particle.coord;
            start_shape.fill_from_normalized(to_cell_coordinates(start),
                                             GHOST_CELLS);
            particle.move(dt);

            Vector3R end = particle.coord;
            end_shape.fill_from_normalized(to_cell_coordinates(end),
                                           start_shape.base_, GHOST_CELLS);
            decompose_esirkepov_current(start_shape, end_shape, qx, qy, qz,
                                        curBuf);

            cellBuf += curBuf;
        }
        int start_x = start_shape.start_.x();
        int start_y = start_shape.start_.y();
        int start_z = start_shape.start_.z();
        flush_current_buffer(fieldJ, cellBuf, start_x, start_y, start_z);
    }
}

void ParticlesArray::correctv_component(const Field3d& fieldE,
                                        const Field3d& fieldEp,
                                        const Field3d& fieldEn,
                                        const Domain& domain, const double dt) {
    if (is_neutral())
        return;

    double jp_cellx = 0;
    double jp_celly = 0;
    double jp_cellz = 0;
    const Bounds bounds = domain.get_bounds();

#pragma omp parallel for reduction(+ : jp_cellx) reduction(+ : jp_celly) reduction(+:jp_cellz)
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto initVelocity = particle.initVelocity;
            const auto velocity = particle.velocity;

            Vector3R end = particle.coord;
            Vector3R coord = end - 0.5 * dt * velocity;

            const Vector3R Ep =
                interpolateE(fieldEp, to_cell_coordinates(coord), SHAPE);
            Vector3R E = interpolateE(fieldE, to_cell_coordinates(coord), SHAPE);
            E += Ep;

            Vector3R v12 = 0.5 * (velocity + initVelocity);
            Vector3R vE =
                Vector3R(v12.x() * E.x(), v12.y() * E.y(), v12.z() * E.z());

            jp_cellx += (0.5 * _mpw * charge) * vE.x();
            jp_celly += (0.5 * _mpw * charge) * vE.y();
            jp_cellz += (0.5 * _mpw * charge) * vE.z();
        }
    }

    const Vector3R energyJeEn =
        calc_JE_component(fieldEn, currentOnGrid, bounds);
    const Vector3R energyJeE =
        calc_JE_component(fieldE, currentOnGrid, bounds);
    const Vector3R energyK = get_kinetic_energy_component();
    Vector3R lambda;
    lambda.x() =
        sqrt(1 + dt * (0.5 * (energyJeEn.x() + energyJeE.x()) - jp_cellx) /
                     energyK.x());
    lambda.y() =
        sqrt(1 + dt * (0.5 * (energyJeEn.y() + energyJeE.y()) - jp_celly) /
                     energyK.y());
    lambda.z() =
        sqrt(1 + dt * (0.5 * (energyJeEn.z() + energyJeE.z()) - jp_cellz) /
                     energyK.z());
    // double lambda2 =
    //     sqrt(1 + Dt *
    //                  (0.5 * (energyJeEn.x() + energyJeE.x() + energyJeEn.y()
    //                  +
    //                          energyJeE.y() + energyJeEn.z() + energyJeE.z())
    //                          -
    //                   jp_cellx - jp_celly - jp_cellz) /
    //                  (energyK.x() + energyK.y() + energyK.z()));

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto velocity = particle.velocity;

            particle.velocity.x() = lambda.x() * velocity.x();
            particle.velocity.y() = lambda.y() * velocity.y();
            particle.velocity.z() = lambda.z() * velocity.z();
            // particle.velocity = lambda2 * velocity;
        }
    }

    // const double energyK2 = get_kinetic_energy();
    std::cout << "lambda " << lambda << "\n";
}

// Very slow function. Fill Lmatrix by each particles
void ParticlesArray::fill_matrixL(Mesh& mesh, const Field3d& fieldB,
                                  const Domain& domain, const double dt,
                                  ShapeType type) {
    if (is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            fill_matrixL_impl_ngp(mesh, fieldB, domain, dt);
            break;
        case ShapeType::Linear:
            fill_matrixL_impl_linear(mesh, fieldB, domain, dt);
            break;
        case ShapeType::Quadratic:
            std::cout << "Fill Lmatrix for quadratic shape function is not implemented" << std::endl;
            exit(-1);    
        
    }
}
// Very slow function. Fill Lmatrix by each particles
void ParticlesArray::fill_matrixL2(Mesh& mesh, const Field3d& fieldB,
                                  const Domain& domain, const double dt,
                                  ShapeType type) {
    if (is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            fill_matrixL_impl_ngp2(mesh, fieldB, domain, dt);
            break;
        case ShapeType::Linear:
            fill_matrixL_impl_linear2(mesh, fieldB, domain, dt);
            break;
        case ShapeType::Quadratic:
            std::cout << "Fill Lmatrix for quadratic shape function is not "
                         "implemented"
                      << std::endl;
            exit(-1);
    }
}
void ParticlesArray::fill_matrixL_impl_ngp(Mesh& mesh, const Field3d& fieldB,
                                           const Domain& domain,
                                           const double dt) {
    constexpr int CHESS_STEP = 3;
#pragma omp parallel
    {
        for (int xStep = 0; xStep < CHESS_STEP; xStep++) {
            for (int yStep = 0; yStep < CHESS_STEP; yStep++) {
#pragma omp for collapse(2) schedule(dynamic, 32)
                for (int ix = xStep; ix < particlesData.size().x();
                     ix += CHESS_STEP) {
                    for (int iy = yStep; iy < particlesData.size().y();
                         iy += CHESS_STEP) {
                        for (int iz = 0; iz < particlesData.size().z(); ++iz) {
                            for (auto& particle : particlesData(ix, iy, iz)) {
                                const auto coord = particle.coord;
                                mesh.update_LmatNGP(coord, domain, charge,
                                                    _mass, _mpw, fieldB, dt);
                            }
                        }
                    }
                }
#pragma omp barrier
            }
        }
    }
}

void ParticlesArray::fill_matrixL_impl_ngp2(Mesh& mesh, const Field3d& fieldB,
                                              const Domain& domain,
                                              const double dt) {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto coord = particle.coord;
            mesh.update_Lmat2_NGP(coord, domain, charge, _mass, _mpw, fieldB, dt);
        }
    }
}

void ParticlesArray::fill_matrixL_impl_linear(Mesh& mesh, const Field3d& fieldB,
                                              const Domain& domain,
                                              const double dt) {
    constexpr int CHESS_STEP = 3;
#pragma omp parallel
    {
        for (int xStep = 0; xStep < CHESS_STEP; xStep++) {
            for (int yStep = 0; yStep < CHESS_STEP; yStep++) {
#pragma omp for collapse(2) schedule(dynamic,32)
                for (int ix = xStep; ix < particlesData.size().x();
                     ix += CHESS_STEP) {
                    for (int iy = yStep; iy < particlesData.size().y();
                         iy += CHESS_STEP) {
                        for (int iz = 0; iz < particlesData.size().z(); ++iz) {
                            for (auto& particle : particlesData(ix, iy, iz)) {
                                const auto coord = particle.coord;
                                mesh.update_Lmat(coord, domain, charge, _mass,
                                                 _mpw, fieldB, dt);
                            }
                        }
                    }
                }
#pragma omp barrier
            }
        }
    }
}

void ParticlesArray::fill_matrixL_impl_linear2(Mesh& mesh, const Field3d& fieldB,
                                              const Domain& domain,
                                              const double dt) {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto coord = particle.coord;
            mesh.update_Lmat2(coord, domain, charge, _mass, _mpw, fieldB, dt);
        }
    }
}
/// @note If shift is false we'll get shape[x - i], shape[x - (i + 0.5)]
/// otherwise
void ParticlesArray::fill_shape(const Node& node, ShapeK& shape,
                                bool isShift) const {
    // Grid coordinates, x,y,z
    double g_x, g_y, g_z;
    const Vector3I g = (isShift) ? node.g05 : node.g;
    const Vector3R r = (isShift) ? node.r - Vector3R(0.5, 0.5, 0.5) : node.r;

#pragma omp simd collapse(3)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            for (int x = 0; x < SHAPE_SIZE; ++x) {
                g_x = g.x() + x;
                g_y = g.y() + y;
                g_z = g.z() + z;

                shape(x, X) = Shape(r.x() - g_x);
                shape(y, Y) = Shape(r.y() - g_y);
                shape(z, Z) = Shape(r.z() - g_z);
            }
        }
    }
}

/// @note If shift is false we'll get shape[x - i], shape[x - (i + 0.5)]
/// otherwise
void ParticlesArray::fill_shape_from_coord(const Vector3R& __r, ShapeK& shape,
                                           bool isShift) const {
    const Vector3R r = to_cell_coordinates(__r);
    const Vector3R rr = (isShift) ? r - Vector3R(0.5, 0.5, 0.5) : r;

    const Vector3I voxel(double_to_int(rr.x()), double_to_int(rr.y()),
                     double_to_int(rr.z()));
    fill_shape(voxel, rr, shape);
}
void ParticlesArray::fill_shape_from_voxel_and_coord(const Vector3I& voxel,
                                                     const Vector3R& __r,
                                                     ShapeK& shape,
                                                     bool isShift) const {
    Vector3R r = to_cell_coordinates(__r);
    const Vector3R rr = (isShift) ? r - Vector3R(0.5, 0.5, 0.5) : r;
    fill_shape(voxel, rr, shape);
}

// voxel coord correspods to particles coord. cell coord is voxel + ghost cells
void ParticlesArray::fill_shape(const Vector3I& voxel, const Vector3R& r,
                                ShapeK& shape) const {
    // Grid coordinates, x,y,z
#pragma omp simd collapse(3)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            for (int x = 0; x < SHAPE_SIZE; ++x) {
                shape(x, X) = Shape(r.x() - (voxel.x() + x));
                shape(y, Y) = Shape(r.y() - (voxel.y() + y));
                shape(z, Z) = Shape(r.z() - (voxel.z() + z));
            }
        }
    }
    shape.cell = voxel + Vector3I(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);
}

Vector3R interpolateE_Chen(const Field3d& fieldE,
                                          const Node& node, ShapeK& sh,
                                          ShapeK& sh_n) {
    Vector3R E = Vector3R(0, 0, 0);

#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            const int ix05 = node.g.x() + GHOST_CELLS;
            const int iy = node.g.y() + y + GHOST_CELLS;
            const int iz = node.g.z() + z + GHOST_CELLS;
            const double w =
                (1. / 3.) * (sh(y, Y) * sh(z, Z) + sh_n(y, Y) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(y, Y) * sh(z, Z) + sh(y, Y) * sh_n(z, Z));
            E.x() += fieldE(ix05, iy, iz, X) * w;
        }
    }
#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            const int ix = node.g.x() + x + GHOST_CELLS;
            const int iy05 = node.g.y() + GHOST_CELLS;
            const int iz = node.g.z() + z + GHOST_CELLS;
            const double w =
                (1. / 3.) * (sh(x, X) * sh(z, Z) + sh_n(x, X) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(x, X) * sh(z, Z) + sh(x, X) * sh_n(z, Z));
            E.y() += fieldE(ix, iy05, iz, Y) * w;
        }
    }
#pragma omp simd collapse(2)
    for (int y = 0; y < SHAPE_SIZE; ++y) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            const int ix = node.g.x() + x + GHOST_CELLS;
            const int iy = node.g.y() + y + GHOST_CELLS;
            const int iz05 = node.g.z() + GHOST_CELLS;
            const double w =
                (1. / 3.) * (sh(x, X) * sh(y, Y) + sh_n(x, X) * sh_n(y, Y)) +
                (1. / 6.) * (sh_n(x, X) * sh(y, Y) + sh(x, X) * sh_n(y, Y));
            E.z() += fieldE(ix, iy, iz05, Z) * w;
        }
    }
    return E;
}

Vector3R interpolateE_Chen(const Field3d& fieldE, ShapeK& sh,
                                          ShapeK& sh_n) {
    Vector3R E = Vector3R(0, 0, 0);

    // check that sh.cell is the same as sh_n.cell
    const Vector3I cell = sh.cell;
#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            const double w =
                (1. / 3.) * (sh(y, Y) * sh(z, Z) + sh_n(y, Y) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(y, Y) * sh(z, Z) + sh(y, Y) * sh_n(z, Z));
            const int ix05 = cell.x();
            const int iy = cell.y() + y;
            const int iz = cell.z() + z;
            E.x() += fieldE(ix05, iy, iz, X) * w;
        }
    }
#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            const double w =
                (1. / 3.) * (sh(x, X) * sh(z, Z) + sh_n(x, X) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(x, X) * sh(z, Z) + sh(x, X) * sh_n(z, Z));
            const int ix = cell.x() + x;
            const int iy05 = cell.y();
            const int iz = cell.z() + z;
            E.y() += fieldE(ix, iy05, iz, Y) * w;
        }
    }
#pragma omp simd collapse(2)
    for (int y = 0; y < SHAPE_SIZE; ++y) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            const double w =
                (1. / 3.) * (sh(x, X) * sh(y, Y) + sh_n(x, X) * sh_n(y, Y)) +
                (1. / 6.) * (sh_n(x, X) * sh(y, Y) + sh(x, X) * sh_n(y, Y));
            const int ix = cell.x() + x;
            const int iy = cell.y() + y;
            const int iz05 = cell.z();
            E.z() += fieldE(ix, iy, iz05, Z) * w;
        }
    }
    return E;
}

Vector3R interpolateE(const Field3d& fieldE, const Node& node,
                                     ShapeK& no, ShapeK& sh) {
    Vector3R E = Vector3R(0, 0, 0);

#pragma omp simd collapse(3)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            for (int x = 0; x < SHAPE_SIZE; ++x) {
                const int ix = node.g.x() + x + GHOST_CELLS;
                const int iy = node.g.y() + y + GHOST_CELLS;
                const int iz = node.g.z() + z + GHOST_CELLS;
                const int ix05 = node.g05.x() + x + GHOST_CELLS;
                const int iy05 = node.g05.y() + y + GHOST_CELLS;
                const int iz05 = node.g05.z() + z + GHOST_CELLS;
                E.x() +=
                    fieldE(ix05, iy, iz, X) * sh(x, X) * no(y, Y) * no(z, Z);
                E.y() +=
                    fieldE(ix, iy05, iz, Y) * no(x, X) * sh(y, Y) * no(z, Z);
                E.z() +=
                    fieldE(ix, iy, iz05, Z) * no(x, X) * no(y, Y) * sh(z, Z);
            }
        }
    }
    return E;
}

Vector3R interpolateB(const Field3d& fieldB, const Node& node,
                                     ShapeK& no, ShapeK& sh) {
    Vector3R B = Vector3R(0, 0, 0);

#pragma omp simd collapse(3)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            for (int x = 0; x < SHAPE_SIZE; ++x) {
                const int ix = node.g.x() + x + GHOST_CELLS;
                const int iy = node.g.y() + y + GHOST_CELLS;
                const int iz = node.g.z() + z + GHOST_CELLS;
                const int ix05 = node.g05.x() + x + GHOST_CELLS;
                const int iy05 = node.g05.y() + y + GHOST_CELLS;
                const int iz05 = node.g05.z() + z + GHOST_CELLS;
                B.x() +=
                    fieldB(ix, iy05, iz05, X) * no(x, X) * sh(y, Y) * sh(z, Z);
                B.y() +=
                    fieldB(ix05, iy, iz05, Y) * sh(x, X) * no(y, Y) * sh(z, Z);
                B.z() +=
                    fieldB(ix05, iy05, iz, Z) * sh(x, X) * sh(y, Y) * no(z, Z);
            }
        }
    }
    return B;
}

// sh for nodes and sh05 for shifted nodes
Vector3R interpolateB(const Field3d& fieldB, ShapeK& sh,
                                     ShapeK& sh05) {
    Vector3R B = Vector3R(0, 0, 0);
#pragma omp simd collapse(3)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            for (int x = 0; x < SHAPE_SIZE; ++x) {
                const int ix = sh.cell.x() + x;
                const int iy = sh.cell.y() + y;
                const int iz = sh.cell.z() + z;
                const int ix05 = sh05.cell.x() + x;
                const int iy05 = sh05.cell.y() + y;
                const int iz05 = sh05.cell.z() + z;
                B.x() += fieldB(ix, iy05, iz05, X) * sh(x, X) * sh05(y, Y) *
                         sh05(z, Z);
                B.y() += fieldB(ix05, iy, iz05, Y) * sh05(x, X) * sh(y, Y) *
                         sh05(z, Z);
                B.z() += fieldB(ix05, iy05, iz, Z) * sh05(x, X) * sh05(y, Y) *
                         sh(z, Z);
            }
        }
    }
    return B;
}

bool ParticlesArray::is_voxel_in_area(const Vector3I& voxel) {
    Vector3I numCells(xCellCount, yCellCount, zCellCount);

    for (int dim = 0; dim < 3; ++dim) {
        if (voxel[dim] < 0 || voxel[dim] >= numCells[dim]) {
            return false;
        }
    }
    return true;
}

bool ParticlesArray::make_periodic_bound_force(Vector3R& point) {
    bool make_bound = false;
    double eps = 1.e-9;
    Vector3R size = Vector3R(static_cast<double>(xCellCount),
                           static_cast<double>(yCellCount),
                           static_cast<double>(zCellCount)) *
                   xCellSize;
    Vector3R new_point = point;
    for (int dim = 0; dim < 3; ++dim) {
        // std::cout << "point " << point(dim) << " " << size(dim) << std::endl;
        if (point[dim] <= eps) {
            new_point[dim] = point[dim] + size[dim] - 2 * eps;
            make_bound = true;
        }
        if (size[dim] <= point[dim] + eps) {
            new_point[dim] = point[dim] - size[dim] + 2 * eps;
            make_bound = true;
        }
    }
    if (!make_bound) {
        std::cout << "fail bound point " << point << std::endl;
    }
    point = new_point;
    return make_bound;
}

// input:
// coord - x_n, velocity - v_n, fieldE - E_{n+1/2}, fieldB - B_{n+1/2}
// output:
// coord - x_{n+1}, velocity - v_{n+1}
double ParticlesArray::track_particle(Vector3R& coord, Vector3R& velocity,
                                      const Field3d& fieldE,
                                      const Field3d& fieldB, double dt,
                                      bool& intersect_bound) {
    const Vector3R cellSize(xCellSize, yCellSize, zCellSize);
    const double bin_size = xCellSize;

    ShapeK shape, new_shape;
    const Vector3R ray_start = coord;
    Vector3R coord_rec = coord;
    Vector3R velocity_rec = velocity;
    Vector3R ray_end = ray_start + velocity * dt;
    Vector3R velocity05;
    Vector3R ray_end_inner;
    bool covergence_cycle = false;
    bool fake_covergence_cycle = false;
    int iter = 0;
    double dt1;
    // may be out if xn is not in area and need change dt
    while (!covergence_cycle) {
        dt1 = dt;
        //    std::cout << "interation " << iter << " " << ray_start << " " <<
        //    ray_end
        //             << "\n";
        std::vector<Vector3I> voxels =
            voxel_traversal(ray_start, ray_end, xCellSize);
        Vector3R current_point = ray_start;
        Vector3R next_point;
        ray_end_inner = ray_end;
        Vector3R E = Vector3R(0, 0, 0);
        Vector3R B = Vector3R(0, 0, 0);
        for (size_t i = 0; i < voxels.size(); i++) {
            if (!is_voxel_in_area(voxels[i])) {
                dt1 = dt * (ray_start - current_point).norm() /
                      (ray_start - ray_end).norm();
                ray_end_inner = current_point;
                //  std::cout << dt1 << " " << voxels[i] << " " << ray_start <<
                //  " " << ray_end << "\n";
                break;
            }

            if (i != voxels.size() - 1) {
                const double t = find_ray_voxel_intersection_parameter(
                    ray_start, ray_end, voxels[i], voxels[i + 1], bin_size);
                next_point = get_point_in_ray(ray_start, ray_end, t);
            } else {
                next_point = ray_end;
            }

            // interpolation and current deposition only in the voxel
            fill_shape_from_voxel_and_coord(voxels[i], current_point, shape,
                                            false);
            fill_shape_from_voxel_and_coord(voxels[i], next_point, new_shape,
                                            false);
            double current_length = (current_point - next_point).norm();
            E += current_length * interpolateE_Chen(fieldE, shape, new_shape);
            // std::cout << interpolateE_Chen(fieldE, shape, new_shape) << "\n";
            // B.z() += current_length * (0.2 * (1 - 0.8 *
            // (0.5*(current_point.y()+next_point.y())- 3.4))); B.x() = B.y() =
            // 0.;
            Vector3R mid_point = 0.5 * (current_point + next_point);
            fill_shape_from_coord(mid_point, shape, false);
            fill_shape_from_coord(mid_point, new_shape, true);

            B += current_length * interpolateB(fieldB, shape, new_shape);
            // std::cout << interpolateB(fieldB, node, shape, new_shape).z() -
            //           (0.2 * (1 - 0.8 * (mid_point.y() - 3.4))) << "\n";
            // std::cout << "track: " << current_point << " " << next_point
            //           << std::endl;
            //  std::cout << E << std::endl;
            //      std::cout << "E B current_length " << E << " " <<  B << " "
            //      << current_length << std::endl;

            current_point = next_point;
        }
        const double length = (ray_end_inner - ray_start).norm();
        if (length > 1.e-16) {
            E /= length;
            B /= length;
        } else {
            fill_shape_from_coord(ray_start, shape, false);
            fill_shape_from_coord(ray_start, new_shape, true);
            E = interpolateE_Chen(fieldE, shape, shape);
            B = interpolateB(fieldB, shape, new_shape);
        }
        //  std::cout << "E B length " << E << " " << B << " " << length
        //           << std::endl;

        const double alpha = 0.5 * dt1 * charge / mass();

        Vector3R ap = velocity + alpha * E;
        velocity05 =
            (ap + alpha * ap.cross(B) + alpha * alpha * ap.dot(B) * B) /
            (1. + pow(alpha * B.norm(), 2));
        // std::cout << "B.norm() " << B.norm() << std::endl;
        // std::cout << "znam " << 1. + pow(alpha * B.norm(), 2) << std::endl;
        // std::cout << "velocity05 " << velocity05 << std::endl;
        Vector3R xn = ray_start + dt1 * velocity05;
        // Vector3R ray_end_inner = ray_start + dt1 * velocity05;

        iter++;
        if (iter > 10) {
            fake_covergence_cycle = true;
          //  fake++;
            //   std::cout << "fake " << fake << std::endl;
            break;
        }
        // std::cout << "iter " << iter << " " << ray_start << " " << ray_end <<
        // " " << xn
        //           << " " << ray_end_inner << std::endl;
        //   std::cout << xn << " " << ray_end_inner << " " << velocity05
        //            << std::endl;
        covergence_cycle = (xn - ray_end_inner).norm() < 1e-10;
        // ray_end = xn;
        ray_end = ray_start + dt * velocity05;
    }

    // std::cout << "interations " << iter << "\n";

    velocity = 2. * velocity05 - velocity;

    intersect_bound = dt1 < dt;

    coord = (intersect_bound) ? ray_end_inner : ray_end;

    if (fake_covergence_cycle) {
        dt1 = track_particle(coord_rec, velocity_rec, fieldE, fieldB, 0.5 * dt,
                             intersect_bound);
        coord = coord_rec;
        velocity = velocity_rec;
    }
    return dt1;
}

void ParticlesArray::updateJ_Chen(const Vector3R value, Field3d& fieldJ,
                                  ShapeK& sh, ShapeK& sh_n) {
// check that sh.cell == sh_n.cell
#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            const int ix05 = sh.cell.x();
            const int iy = sh.cell.y() + y;
            const int iz = sh.cell.z() + z;
            const double w =
                (1. / 3.) * (sh(y, Y) * sh(z, Z) + sh_n(y, Y) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(y, Y) * sh(z, Z) + sh(y, Y) * sh_n(z, Z));
            fieldJ(ix05, iy, iz, X) += value.x() * w;
        }
    }
#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            const int ix = sh.cell.x() + x;
            const int iy05 = sh.cell.y();
            const int iz = sh.cell.z() + z;
            const double w =
                (1. / 3.) * (sh(x, X) * sh(z, Z) + sh_n(x, X) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(x, X) * sh(z, Z) + sh(x, X) * sh_n(z, Z));
            fieldJ(ix, iy05, iz, Y) += value.y() * w;
        }
    }
#pragma omp simd collapse(2)
    for (int y = 0; y < SHAPE_SIZE; ++y) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            const int ix = sh.cell.x() + x;
            const int iy = sh.cell.y() + y;
            const int iz05 = sh.cell.z();
            const double w =
                (1. / 3.) * (sh(x, X) * sh(y, Y) + sh_n(x, X) * sh_n(y, Y)) +
                (1. / 6.) * (sh_n(x, X) * sh(y, Y) + sh(x, X) * sh_n(y, Y));
            fieldJ(ix, iy, iz05, Z) += value.z() * w;
        }
    }
}

// input: coord_start - x_n,  coord_end - x_{n+1}, fieldJ -  J_{n+1/2}
void ParticlesArray::calc_current_Chen(const Vector3R& coord_start,
                                       const Vector3R& coord_end,
                                       Field3d& fieldJ, const double dt) {
    const Vector3R cellSize(xCellSize, yCellSize, zCellSize);
    // todo : check xCellSize = yCellSize = zCellSize
    const double bin_size = xCellSize;
    const double volume = 1;   // bin_size * bin_size * bin_size;
    // std::cout << "volume " << volume << std::endl;

    ShapeK shape, new_shape;

    std::vector<Vector3I> voxels =
        voxel_traversal(coord_start, coord_end, xCellSize);
    Vector3R current_point = coord_start;
    Vector3R next_point;

    for (size_t i = 0; i < voxels.size(); i++) {
        // if (!is_voxel_in_area(voxels[i])) {
        //     std::cout << "voxel out of area " << voxels[i] << " " <<
        //     coord_start << " " << coord_end << "\n"; break;
        // }

        if (i != voxels.size() - 1) {
            double t = find_ray_voxel_intersection_parameter(
                coord_start, coord_end, voxels[i], voxels[i + 1], bin_size);
            next_point = get_point_in_ray(coord_start, coord_end, t);
        } else {
            next_point = coord_end;
        }
        // std::cout << current_point << " " << next_point << "\n";
        // Node node(current_point, cellSize);
        // Node new_node(next_point, cellSize);
        fill_shape_from_voxel_and_coord(voxels[i], current_point, shape, false);
        fill_shape_from_voxel_and_coord(voxels[i], next_point, new_shape,
                                        false);
        Vector3R current_length = next_point - current_point;
        Vector3R value = _mpw * charge * (current_length / volume) / dt;
        // std::cout << "cur len "<< current_length << "value " << value <<
        // "\n";

        updateJ_Chen(value, fieldJ, shape, new_shape);

        current_point = next_point;
    }
}

void ParticlesArray::push_Chen(const Field3d& fieldE,
                               const Field3d& fieldB, double dt) {
    if (is_neutral())
        return;

    // std::cout << "fieldB05 " << fieldB05(2,3,4,0)  << "\n";
    currentOnGrid.setZero();

    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            particle.coord = particle.initCoord;
            particle.velocity = particle.initVelocity;
        }
    }
    // const double dt_eps = 1.e-6;
    //  Vector3R cellSize(xCellSize, yCellSize, zCellSize);
    //  double bin_size = xCellSize;

#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < size(); ++pk) {
        // ShapeK shape, new_shape;

        for (auto& particle : particlesData(pk)) {
            // std::cout << counter++ << "\n";
            // std::cout << particle << "\n";
            double dtt = 0;
            double dt2 = dt;
            while (dtt < dt) {
                Vector3R coord = particle.coord;
                Vector3R velocity = particle.velocity;
                bool intersect_bound = false;
                // std::cout << "track start "
                //           << "\n";

                double dt_fact = track_particle(coord, velocity, fieldE, fieldB,
                                                dt2, intersect_bound);
                //      std::cout << "dt_fact " << dt_fact << "\n";
                calc_current_Chen(particle.coord, coord, currentOnGrid, dt);
                //    std::cout << "calc_current_Chen " << "\n";

                // std::cout << "track done " << dt_fact << "\n";
                if (intersect_bound) {
                    //  std::cout << "intersect bound" << particle << std::endl;
                    if (!make_periodic_bound_force(coord)) // todo - if periodic boundary
                        std::cout << particle.coord << " " << dt_fact
                                  << " error in bound correction" << std::endl;
                }
                dtt += dt_fact;

                dt2 = dt - dtt;
                // if (intersect_bound) {std::cout<< dtt << " " << dt2 << "\n";
                // exit(0);}

                particle.coord = coord;
                particle.velocity = velocity;
                //    std::cout << particle << "\n";
            }
        }
    }

}
