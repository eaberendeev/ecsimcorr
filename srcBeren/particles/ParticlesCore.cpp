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

    const double qx = charge * domain_.cell_size().x() / (6 * dt) * mpw_;
    const double qy = charge * domain_.cell_size().y() / (6 * dt) * mpw_;
    const double qz = charge * domain_.cell_size().z() / (6 * dt) * mpw_;

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
        curBuf.zero();
        start_shape.fill_zero();

        for (auto& particle : particles) {
            Vector3R start = particle.coord;
            start_shape.fill_from_normalized(domain_.to_cell_coordinates(start),
                                             GHOST_CELLS);
            particle.move(dt);

            Vector3R end = particle.coord;
            end_shape.fill_from_normalized(domain_.to_cell_coordinates(end),
                                           start_shape.base_, GHOST_CELLS);
            decompose_esirkepov_current(start_shape, end_shape, qx, qy, qz,
                                        curBuf);

            cellBuf += curBuf;
        }
        auto [start_x, start_y, start_z] = start_shape.start_.split();
        flush_current_buffer(fieldJ, cellBuf, start_x, start_y, start_z);
    }
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
                                                    mass_, mpw_, fieldB, dt);
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
            mesh.update_Lmat2_NGP(coord, domain, charge, mass_, mpw_, fieldB, dt);
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
                                mesh.update_Lmat(coord, domain, charge, mass_,
                                                 mpw_, fieldB, dt);
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
            mesh.update_Lmat2(coord, domain, charge, mass_, mpw_, fieldB, dt);
        }
    }
}
