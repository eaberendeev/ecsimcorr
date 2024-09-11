#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "containers.h"
#include "service.h"
#include "voxel_traversal.h"

void ParticlesArray::move(double dt) {
#pragma omp parallel for
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.move(dt);
        }
    }
}

void ParticlesArray::predict_velocity(const Mesh& mesh, const Domain& domain,
                                      const double dt) {
#pragma omp parallel for
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            const auto coord = particle.coord;
            const auto velocity = particle.velocity;
            double3 E = get_fieldE_in_pos(mesh.fieldE, coord, domain);
            const double3 B = get_fieldB_in_pos(mesh.fieldB, coord, domain);
            const double3 En = get_fieldE_in_pos(mesh.fieldEp, coord, domain);
            E = 0.5 * (E + En);
            const auto beta = dt * charge / _mass;
            const auto alpha = 0.5 * beta * mag(B);
            const auto alpha2 = alpha * alpha;
            const auto h = unit(B);

            const auto v12 =
                (1. / (1. + alpha2)) *
                (velocity + alpha * cross(velocity, h) +
                 alpha2 * dot(h, velocity) * h +
                 0.5 * beta *
                     (E + alpha * cross(E, h) + alpha2 * dot(E, h) * h));

            particle.velocity = 2. * v12 - velocity;
        }
    }
}

// TO DO:: merge with void ParticlesArray::move_and_calc_current(const double
// dt, Field3d& fieldJ)
void ParticlesArray::move_and_calc_current(const double dt) {
    constexpr auto SMAX = 2 * SHAPE_SIZE;

    const double conx = xCellSize / (6 * dt) * _mpw;
    const double cony = yCellSize / (6 * dt) * _mpw;
    const double conz = zCellSize / (6 * dt) * _mpw;

#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        for (auto& particle : particlesData(pk)) {
            double3 start = particle.coord;

            particle.move(dt);

            double3 end = particle.coord;

            double xx = start.x() / xCellSize;
            double yy = start.y() / yCellSize;
            double zz = start.z() / zCellSize;

            double xn = end.x() / xCellSize;
            double yn = end.y() / yCellSize;
            double zn = end.z() / zCellSize;

            int xk = int(xx);
            int yk = int(yy);
            int zk = int(zz);

            for (int n = 0; n < SMAX; ++n) {
                for (int m = 0; m < SMAX; ++m) {
                    for (int k = 0; k < SMAX; ++k) {
                        jx[n][m][k] = 0.;
                        jy[n][m][k] = 0.;
                        jz[n][m][k] = 0.;
                    }
                }
            }
            double arg;
            for (int n = 0; n < SMAX; ++n) {
                arg = -xx + double(xk - GHOST_CELLS + n);
                sx[n] = Shape2(arg);
                arg = -yy + double(yk - GHOST_CELLS + n);
                sy[n] = Shape2(arg);
                arg = -zz + double(zk - GHOST_CELLS + n);
                sz[n] = Shape2(arg);
                arg = -xn + double(xk - GHOST_CELLS + n);
                sx_n[n] = Shape2(arg);
                arg = -yn + double(yk - GHOST_CELLS + n);
                sy_n[n] = Shape2(arg);
                arg = -zn + double(zk - GHOST_CELLS + n);
                sz_n[n] = Shape2(arg);
            }

            for (int n = 0; n < SMAX; ++n) {
                int indx = xk + n;
                for (int m = 0; m < SMAX; ++m) {
                    int indy = yk + m;
                    for (int k = 0; k < SMAX; ++k) {
                        int indz = zk + k;

                        if (n == 0)
                            jx[n][m][k] = -charge * conx * (sx_n[n] - sx[n]) *
                                          (sy_n[m] * (2 * sz_n[k] + sz[k]) +
                                           sy[m] * (2 * sz[k] + sz_n[k]));

                        if (n > 0 && n < SMAX - 1)
                            jx[n][m][k] = jx[n - 1][m][k] -
                                          charge * conx * (sx_n[n] - sx[n]) *
                                              (sy_n[m] * (2 * sz_n[k] + sz[k]) +
                                               sy[m] * (2 * sz[k] + sz_n[k]));

                        if (m == 0)
                            jy[n][m][k] = -charge * cony * (sy_n[m] - sy[m]) *
                                          (sx_n[n] * (2 * sz_n[k] + sz[k]) +
                                           sx[n] * (2 * sz[k] + sz_n[k]));
                        if (m > 0 && m < SMAX - 1)
                            jy[n][m][k] = jy[n][m - 1][k] -
                                          charge * cony * (sy_n[m] - sy[m]) *
                                              (sx_n[n] * (2 * sz_n[k] + sz[k]) +
                                               sx[n] * (2 * sz[k] + sz_n[k]));

                        if (k == 0)
                            jz[n][m][k] = -charge * conz * (sz_n[k] - sz[k]) *
                                          (sy_n[m] * (2 * sx_n[n] + sx[n]) +
                                           sy[m] * (2 * sx[n] + sx_n[n]));
                        if (k > 0 && k < SMAX - 1)
                            jz[n][m][k] = jz[n][m][k - 1] -
                                          charge * conz * (sz_n[k] - sz[k]) *
                                              (sy_n[m] * (2 * sx_n[n] + sx[n]) +
                                               sy[m] * (2 * sx[n] + sx_n[n]));

#pragma omp atomic update
                        currentOnGrid(indx, indy, indz, 0) += jx[n][m][k];
#pragma omp atomic update
                        currentOnGrid(indx, indy, indz, 1) += jy[n][m][k];
#pragma omp atomic update
                        currentOnGrid(indx, indy, indz, 2) += jz[n][m][k];
                    }
                }
            }
        }
    }
}

void ParticlesArray::move_and_calc_current(const double dt, Field3d& fieldJ) {
    constexpr auto SMAX = 2 * SHAPE_SIZE;

    const double conx = xCellSize / (6 * dt) * _mpw;
    const double cony = yCellSize / (6 * dt) * _mpw;
    const double conz = zCellSize / (6 * dt) * _mpw;

#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        double arg;
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        for (auto& particle : particlesData(pk)) {
            double3 start = particle.coord;

            particle.move(dt);

            double3 end = particle.coord;

            double xx = start.x() / xCellSize;
            double yy = start.y() / yCellSize;
            double zz = start.z() / zCellSize;

            double xn = end.x() / xCellSize;
            double yn = end.y() / yCellSize;
            double zn = end.z() / zCellSize;

            int xk = int(xx);
            int yk = int(yy);
            int zk = int(zz);

            for (int n = 0; n < SMAX; ++n) {
                for (int m = 0; m < SMAX; ++m) {
                    for (int k = 0; k < SMAX; ++k) {
                        jx[n][m][k] = 0.;
                        jy[n][m][k] = 0.;
                        jz[n][m][k] = 0.;
                    }
                }
            }

            for (int n = 0; n < SMAX; ++n) {
                arg = -xx + double(xk - GHOST_CELLS + n);
                sx[n] = Shape2(arg);
                arg = -yy + double(yk - GHOST_CELLS + n);
                sy[n] = Shape2(arg);
                arg = -zz + double(zk - GHOST_CELLS + n);
                sz[n] = Shape2(arg);
                arg = -xn + double(xk - GHOST_CELLS + n);
                sx_n[n] = Shape2(arg);
                arg = -yn + double(yk - GHOST_CELLS + n);
                sy_n[n] = Shape2(arg);
                arg = -zn + double(zk - GHOST_CELLS + n);
                sz_n[n] = Shape2(arg);
            }

            for (int n = 0; n < SMAX; ++n) {
                int indx = xk + n;
                for (int m = 0; m < SMAX; ++m) {
                    int indy = yk + m;
                    for (int k = 0; k < SMAX; ++k) {
                        int indz = zk + k;

                        if (n == 0)
                            jx[n][m][k] = -charge * conx * (sx_n[n] - sx[n]) *
                                          (sy_n[m] * (2 * sz_n[k] + sz[k]) +
                                           sy[m] * (2 * sz[k] + sz_n[k]));

                        if (n > 0 && n < SMAX - 1)
                            jx[n][m][k] = jx[n - 1][m][k] -
                                          charge * conx * (sx_n[n] - sx[n]) *
                                              (sy_n[m] * (2 * sz_n[k] + sz[k]) +
                                               sy[m] * (2 * sz[k] + sz_n[k]));

                        if (m == 0)
                            jy[n][m][k] = -charge * cony * (sy_n[m] - sy[m]) *
                                          (sx_n[n] * (2 * sz_n[k] + sz[k]) +
                                           sx[n] * (2 * sz[k] + sz_n[k]));
                        if (m > 0 && m < SMAX - 1)
                            jy[n][m][k] = jy[n][m - 1][k] -
                                          charge * cony * (sy_n[m] - sy[m]) *
                                              (sx_n[n] * (2 * sz_n[k] + sz[k]) +
                                               sx[n] * (2 * sz[k] + sz_n[k]));

                        if (k == 0)
                            jz[n][m][k] = -charge * conz * (sz_n[k] - sz[k]) *
                                          (sy_n[m] * (2 * sx_n[n] + sx[n]) +
                                           sy[m] * (2 * sx[n] + sx_n[n]));
                        if (k > 0 && k < SMAX - 1)
                            jz[n][m][k] = jz[n][m][k - 1] -
                                          charge * conz * (sz_n[k] - sz[k]) *
                                              (sy_n[m] * (2 * sx_n[n] + sx[n]) +
                                               sy[m] * (2 * sx[n] + sx_n[n]));

#pragma omp atomic update
                        fieldJ(indx, indy, indz, 0) += jx[n][m][k];
#pragma omp atomic update
                        fieldJ(indx, indy, indz, 1) += jy[n][m][k];
#pragma omp atomic update
                        fieldJ(indx, indy, indz, 2) += jz[n][m][k];
                    }
                }
            }
        }
    }
}

void ParticlesArray::calc_Esirkepov_current(const double dt,
                                            Field3d& fieldJ) const {
    constexpr auto SMAX = 2 * SHAPE_SIZE;
    const double conx = xCellSize / (6 * dt) * _mpw;
    const double cony = yCellSize / (6 * dt) * _mpw;
    const double conz = zCellSize / (6 * dt) * _mpw;
#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        double arg;
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        double q = charge;

        for (auto& particle : particlesData(pk)) {
            const auto start = particle.initCoord;
            const auto end = particle.coord;

            const auto xx = start.x() / xCellSize;
            const auto yy = start.y() / yCellSize;
            const auto zz = start.z() / zCellSize;

            const auto xn = end.x() / xCellSize;
            const auto yn = end.y() / yCellSize;
            const auto zn = end.z() / zCellSize;

            const auto xk = int(xx);
            const auto yk = int(yy);
            const auto zk = int(zz);

            for (int n = 0; n < SMAX; ++n) {
                for (int m = 0; m < SMAX; ++m) {
                    for (int k = 0; k < SMAX; ++k) {
                        jx[n][m][k] = 0.;
                        jy[n][m][k] = 0.;
                        jz[n][m][k] = 0.;
                    }
                }
            }

            for (int n = 0; n < SMAX; ++n) {
                arg = -xx + double(xk - GHOST_CELLS + n);
                sx[n] = Shape(arg);
                arg = -yy + double(yk - GHOST_CELLS + n);
                sy[n] = Shape(arg);
                arg = -zz + double(zk - GHOST_CELLS + n);
                sz[n] = Shape(arg);
                arg = -xn + double(xk - GHOST_CELLS + n);
                sx_n[n] = Shape(arg);
                arg = -yn + double(yk - GHOST_CELLS + n);
                sy_n[n] = Shape(arg);
                arg = -zn + double(zk - GHOST_CELLS + n);
                sz_n[n] = Shape(arg);
            }

            for (int n = 0; n < SMAX; ++n) {
                const auto indx = xk + n;
                for (int m = 0; m < SMAX; ++m) {
                    const auto indy = yk + m;
                    for (int k = 0; k < SMAX; ++k) {
                        const auto indz = zk + k;

                        if (n == 0)
                            jx[n][m][k] = -q * conx * (sx_n[n] - sx[n]) *
                                          (sy_n[m] * (2 * sz_n[k] + sz[k]) +
                                           sy[m] * (2 * sz[k] + sz_n[k]));

                        if (n > 0 && n < SMAX - 1)
                            jx[n][m][k] = jx[n - 1][m][k] -
                                          q * conx * (sx_n[n] - sx[n]) *
                                              (sy_n[m] * (2 * sz_n[k] + sz[k]) +
                                               sy[m] * (2 * sz[k] + sz_n[k]));

                        if (m == 0)
                            jy[n][m][k] = -q * cony * (sy_n[m] - sy[m]) *
                                          (sx_n[n] * (2 * sz_n[k] + sz[k]) +
                                           sx[n] * (2 * sz[k] + sz_n[k]));
                        if (m > 0 && m < SMAX - 1)
                            jy[n][m][k] = jy[n][m - 1][k] -
                                          q * cony * (sy_n[m] - sy[m]) *
                                              (sx_n[n] * (2 * sz_n[k] + sz[k]) +
                                               sx[n] * (2 * sz[k] + sz_n[k]));

                        if (k == 0)
                            jz[n][m][k] = -q * conz * (sz_n[k] - sz[k]) *
                                          (sy_n[m] * (2 * sx_n[n] + sx[n]) +
                                           sy[m] * (2 * sx[n] + sx_n[n]));
                        if (k > 0 && k < SMAX - 1)
                            jz[n][m][k] = jz[n][m][k - 1] -
                                          q * conz * (sz_n[k] - sz[k]) *
                                              (sy_n[m] * (2 * sx_n[n] + sx[n]) +
                                               sy[m] * (2 * sx[n] + sx_n[n]));
#pragma omp atomic update
                        fieldJ(indx, indy, indz, 0) += jx[n][m][k];
#pragma omp atomic update
                        fieldJ(indx, indy, indz, 1) += jy[n][m][k];
#pragma omp atomic update
                        fieldJ(indx, indy, indz, 2) += jz[n][m][k];
                    }
                }
            }
        }
    }
}

void ParticlesArray::correctv(Mesh& mesh, const Domain& domain,
                              const double dt) {
    std::array<double, 20> ldistr;
    for (auto& val : ldistr) {
        val = 0.0;
    }

    double jp_cell = 0;
#pragma omp parallel for reduction(+ : jp_cell)
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto initVelocity = particle.initVelocity;
            const auto velocity = particle.velocity;

            double3 end = particle.coord;
            double3 coord = end - 0.5 * dt * velocity;
            double3 Ep = get_fieldE_in_pos(mesh.fieldEp, coord, domain);
            double3 E = get_fieldE_in_pos(mesh.fieldE, coord, domain);

            double3 v12 = 0.5 * (velocity + initVelocity);

            jp_cell += 0.5 * _mpw * charge * dot(v12, (Ep + E));
        }
    }

    mesh.fieldJp_full.data() =
        mesh.fieldJp.data() +
        mesh.Lmat2 * (mesh.fieldE.data() + mesh.fieldEp.data()) / dt;

    const double energyJeEn = mesh.calc_JE(mesh.fieldEn, currentOnGrid);
    const double energyJeE = mesh.calc_JE(mesh.fieldE, currentOnGrid);
    const double energyJpEp = mesh.calc_JE(mesh.fieldEp, mesh.fieldJp_full);
    const double energyJpE = mesh.calc_JE(mesh.fieldE, mesh.fieldJp_full);
    const double energyK = get_kinetic_energy();
    const double lambda =
        sqrt(1 + dt * (0.5 * (energyJeEn + energyJeE) - jp_cell) / energyK);

#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto velocity = particle.velocity;

            particle.velocity = lambda * velocity;
        }
    }

    const double energyK2 = get_kinetic_energy();
    std::cout << "lambda " << lambda << " " << lambda * lambda << " "
              << energyK2 - energyK << " "
              << 0.5 * dt * (energyJeEn + energyJeE - energyJpEp - energyJpE)
              << "\n";
}

void ParticlesArray::correctv_component(Mesh& mesh, const Domain& domain,
                                        const double dt) {
    double jp_cellx = 0;
    double jp_celly = 0;
    double jp_cellz = 0;
#pragma omp parallel for reduction(+ : jp_cellx) reduction(+ : jp_celly) reduction(+:jp_cellz)
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto initVelocity = particle.initVelocity;
            const auto velocity = particle.velocity;

            double3 end = particle.coord;
            double3 coord = end - 0.5 * dt * velocity;

            double3 Ep = get_fieldE_in_pos(mesh.fieldEp, coord, domain);
            double3 E = get_fieldE_in_pos(mesh.fieldE, coord, domain);
            E += Ep;

            double3 v12 = 0.5 * (velocity + initVelocity);
            double3 vE =
                double3(v12.x() * E.x(), v12.y() * E.y(), v12.z() * E.z());

            jp_cellx += (0.5 * _mpw * charge) * vE.x();
            jp_celly += (0.5 * _mpw * charge) * vE.y();
            jp_cellz += (0.5 * _mpw * charge) * vE.z();
        }
    }

    const double3 energyJeEn =
        mesh.calc_JE_component(mesh.fieldEn, currentOnGrid);
    const double3 energyJeE =
        mesh.calc_JE_component(mesh.fieldE, currentOnGrid);
    const double3 energyK = get_kinetic_energy_component();
    double3 lambda;
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

#pragma omp parallel for
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

void ParticlesArray::predict_current(const Field3d& fieldB, Field3d& fieldJ,
                                     const Domain& domain, const double dt) {
    constexpr auto SMAX = SHAPE_SIZE;
#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        int i, j, k;
        int i05, j05, k05;
        double qp = charge;

        alignas(64) double wx[SMAX], wy[SMAX], wz[SMAX];
        alignas(64) double wx05[SMAX], wy05[SMAX], wz05[SMAX];

        for (auto& particle : particlesData(pk)) {
            double3 coord = particle.coord;
            double3 velocity = particle.velocity;

            double x = coord.x() / xCellSize + GHOST_CELLS;
            double y = coord.y() / yCellSize + GHOST_CELLS;
            double z = coord.z() / zCellSize + GHOST_CELLS;
            double x05 = x - 0.5;
            double y05 = y - 0.5;
            double z05 = z - 0.5;

            const auto ix = int(x);
            const auto iy = int(y);
            const auto iz = int(z);
            const auto ix05 = int(x05);
            const auto iy05 = int(y05);
            const auto iz05 = int(z05);

            wx[1] = (x - ix);
            wx[0] = 1 - wx[1];
            wy[1] = (y - iy);
            wy[0] = 1 - wy[1];
            wz[1] = (z - iz);
            wz[0] = 1 - wz[1];

            wx05[1] = (x05 - ix05);
            wx05[0] = 1 - wx05[1];
            wy05[1] = (y05 - iy05);
            wy05[0] = 1 - wy05[1];
            wz05[1] = (z05 - iz05);
            wz05[0] = 1 - wz05[1];

            double3 B = get_fieldB_in_pos(fieldB, coord, domain);

            double beta = dt * qp / _mass;
            double alpha = 0.5 * beta * mag(B);
            double alpha2 = alpha * alpha;
            double3 h = unit(B);

            double3 current = qp * _mpw / (1. + alpha2) *
                              (velocity + alpha * cross(velocity, h) +
                               alpha2 * dot(h, velocity) * h);

            for (int nx = 0; nx < SMAX; ++nx) {
                i = ix + nx;
                i05 = ix05 + nx;
                for (int ny = 0; ny < SMAX; ++ny) {
                    j = iy + ny;
                    j05 = iy05 + ny;
                    for (int nz = 0; nz < SMAX; ++nz) {
                        k = iz + nz;
                        k05 = iz05 + nz;
                        double sx = wx05[nx] * wy[ny] * wz[nz];
                        double sy = wx[nx] * wy05[ny] * wz[nz];
                        double sz = wx[nx] * wy[ny] * wz05[nz];
#pragma omp atomic update
                        fieldJ(i05, j, k, 0) += sx * current.x();
#pragma omp atomic update
                        fieldJ(i, j05, k, 1) += sy * current.y();
#pragma omp atomic update
                        fieldJ(i, j, k05, 2) += sz * current.z();
                    }
                }
            }
        }
    }
}

// Very slow function. Fill Lmatrix by each particles
void ParticlesArray::get_L(Mesh& mesh, const Domain& domain, const double dt) {
#pragma omp parallel
    {
        for (int xStep = 0; xStep < 4; xStep++) {
            for (int yStep = 0; yStep < 4; yStep++) {
#pragma omp for collapse(2) schedule(dynamic)
                for (int ix = xStep; ix < particlesData.size().x(); ix += 4) {
                    for (int iy = yStep; iy < particlesData.size().y();
                         iy += 4) {
                        for (int iz = 0; iz < particlesData.size().z(); ++iz) {
                            for (auto& particle : particlesData(ix, iy, iz)) {
                                const auto coord = particle.coord;
                                mesh.update_Lmat(coord, domain, charge, _mass,
                                                 _mpw, dt);
                            }
                        }
                    }
                }
#pragma omp barrier
            }
        }
    }
}


/// @note If shift is false we'll get shape[x - i], shape[x - (i + 0.5)]
/// otherwise
void ParticlesArray::fill_shape(const Node& node, ShapeK& shape,
                                bool isShift) const {
    // Grid coordinates, x,y,z
    double g_x, g_y, g_z;
    const int3 g = (isShift) ? node.g05 : node.g;
    const double3 r = (isShift) ? node.r - double3(0.5, 0.5, 0.5) : node.r;

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

double3 ParticlesArray::interpolateE_Chen(const Field3d& fieldE,
                                          const Node& node, ShapeK& sh,
                                          ShapeK& sh_n) {
    double3 E = double3(0, 0, 0);

#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            int ix05 = node.g.x() + GHOST_CELLS;
            int iy = node.g.y() + y + GHOST_CELLS;
            int iz = node.g.z() + z + GHOST_CELLS;
            double w =
                (1. / 3.) * (sh(y, Y) * sh(z, Z) + sh_n(y, Y) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(y, Y) * sh(z, Z) + sh(y, Y) * sh_n(z, Z));
            E.x() += fieldE(ix05, iy, iz, X) * w;
        }
    }
#pragma omp simd collapse(2)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            int ix = node.g.x() + x + GHOST_CELLS;
            int iy05 = node.g.y() + GHOST_CELLS;
            int iz = node.g.z() + z + GHOST_CELLS;
            double w =
                (1. / 3.) * (sh(x, X) * sh(z, Z) + sh_n(x, X) * sh_n(z, Z)) +
                (1. / 6.) * (sh_n(x, X) * sh(z, Z) + sh(x, X) * sh_n(z, Z));
            E.y() += fieldE(ix, iy05, iz, Y) * w;
        }
    }
#pragma omp simd collapse(2)
    for (int y = 0; y < SHAPE_SIZE; ++y) {
        for (int x = 0; x < SHAPE_SIZE; ++x) {
            int ix = node.g.x() + x + GHOST_CELLS;
            int iy = node.g.y() + y + GHOST_CELLS;
            int iz05 = node.g.z() + GHOST_CELLS;
            double w =
                (1. / 3.) * (sh(x, X) * sh(y, Y) + sh_n(x, X) * sh_n(y, Y)) +
                (1. / 6.) * (sh_n(x, X) * sh(y, Y) + sh(x, X) * sh_n(y, Y));
            E.z() += fieldE(ix, iy, iz05, Z) * w;
        }
    }
    return E;
}

double3 ParticlesArray::interpolateE(const Field3d& fieldE, const Node& node,
                                     ShapeK& no, ShapeK& sh) {
    double3 E = double3(0, 0, 0);

#pragma omp simd collapse(3)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            for (int x = 0; x < SHAPE_SIZE; ++x) {
                int ix = node.g.x() + x + GHOST_CELLS;
                int iy = node.g.y() + y + GHOST_CELLS;
                int iz = node.g.z() + z + GHOST_CELLS;
                int ix05 = node.g05.x() + x + GHOST_CELLS;
                int iy05 = node.g05.y() + y + GHOST_CELLS;
                int iz05 = node.g05.z() + z + GHOST_CELLS;
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

double3 ParticlesArray::interpolateB(const Field3d& fieldB, const Node& node,
                                     ShapeK& no, ShapeK& sh) {
    double3 B = double3(0, 0, 0);

#pragma omp simd collapse(3)
    for (int z = 0; z < SHAPE_SIZE; ++z) {
        for (int y = 0; y < SHAPE_SIZE; ++y) {
            for (int x = 0; x < SHAPE_SIZE; ++x) {
                int ix = node.g.x() + x + GHOST_CELLS;
                int iy = node.g.y() + y + GHOST_CELLS;
                int iz = node.g.z() + z + GHOST_CELLS;
                int ix05 = node.g05.x() + x + GHOST_CELLS;
                int iy05 = node.g05.y() + y + GHOST_CELLS;
                int iz05 = node.g05.z() + z + GHOST_CELLS;
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

// bool ParticlesArray::boundary_correction(double3& coord) {
//     for (int i = 0; i < 3; ++i) {
//         if (boundary_correction(coord(i), length, i)) {
//             return true;
//         }
//     }
//     return false;
// }

// bool ParticlesArray::boundary_correction(double& coord, const int dim){
//     const double length = cellSize(dim) * cellCount(dim);
//     if (coord(dim) < 0) {
//         coord(dim) += length;
//         return true;
//     }
//     if (coord(dim) > length) {
//         coord(dim) -= length;
//         return true;
//     }
//     return false;
// }

// void ParticlesArray::make_periodic_bound(int3& voxel) {
//     int3 cellCount(xCellCount, yCellCount, zCellCount);
//     for (int dim = 0; dim < 3; ++dim) {
//         if (voxel(dim) < 0) {
//             voxel(dim) += cellCount(dim);
//         }
//         if (voxel(dim) >= cellCount(dim)) {
//             voxel(dim) -= cellCount(dim);
//         }
//     }
// }

// bool ParticlesArray::make_periodic_bound(int3& voxel, double3& ray_start,
//                                          double3& ray_end) {
//     bool is_done = true;
//     int3 cellCount(xCellCount, yCellCount, zCellCount);
//     double bin_size = xCellSize;
//     for (int dim = 0; dim < 3; ++dim) {
//         if (voxel(dim) < 0) {
//             voxel(dim) += cellCount(dim);
//             ray_start(dim) += cellCount(dim) * bin_size;
//             ray_end(dim) += cellCount(dim) * bin_size;
//             return is_done;
//         }
//         if (voxel(dim) >= cellCount(dim)) {
//             voxel(dim) -= cellCount(dim);
//             ray_start(dim) -= cellCount(dim) * bin_size;
//             ray_end(dim) -= cellCount(dim) * bin_size;
//             return is_done;
//         }
//     }
//     return !is_done;
// }

// bool ParticlesArray::in_exended_domain(const double3& coord) {
//     double3 upper_bound =
//         double3(xCellCount + 0.5, yCellCount + 0.5, zCellCount + 0.5) *
//         xCellSize;
//     double3 lower_bound = -double3(0.5, 0.5, 0.5) * xCellSize;
//     for (int i = 0; i < 3; ++i) {
//         if (coord(i) < lower_bound(i) || coord(i) > upper_bound(i)) {
//             return false;
//         }
//     }
//     return true;
// }
// void track_particle(const double3& point_start, const double3& point_end, double3& E, double3& B, const double length) {
//     bool is_move = true;
//     double3 ray_start = point_start;
//     double3 ray_end = point_end;
//     double bin_size = xCellSize;

//     while (is_move) {
//         std::vector<int3> voxels =
//             voxel_traversal(ray_start, ray_end, bin_size);
//         double3 current_point = ray_start;
//         double3 next_point;
//         int3 final_voxel = voxels[voxels.size() - 1];
//         make_periodic_bound(final_voxel);
//         for (int i = 0; i < voxels.size(); i++) {
//             if (i != voxels.size() - 1) {
//                 double t = find_ray_voxel_intersection_parameter(
//                     ray_start, ray_end, voxels[i], voxels[i + 1], bin_size);
//                 next_point = get_point_in_ray(ray_start, ray_end, t);
//             } else {
//                 next_point = ray_end;
//             }
//             Node node(current_point, bin_size);
//             Node new_node(next_point, bin_size);
//             fill_shape(node, shape, false);
//             fill_shape(new_node, new_shape, false);
//             double current_length = (current_point - next_point).norm();
//             E += (current_length / length) *
//                  interpolateE_Chen(fieldE, node, shape, new_shape);
//             B += (current_length / length) *
//                  interpolateB(fieldB, node, shape, new_shape);
//             std::cout << "track: " << current_point << " " << next_point
//                       << std::endl;
//             // std::cout << E << std::endl;
//             current_point = next_point;
//             if (make_periodic_bound(voxels[i + 1], current_point, ray_end)) {
//                 ray_start = current_point;
//                 break;
//             }
//             is_move = (voxels[i] != final_voxel);
//         }
//     }
// }

// void ParticlesArray::push_Chen(const Field3d& fieldE, const Field3d& fieldB,
//                                double dt) {
//     double3 cellSize(xCellSize, yCellSize, zCellSize);
//     double bin_size = xCellSize;
// #pragma omp parallel for
//     for (auto pk = 0; pk < size(); ++pk) {
//         ShapeK shape, new_shape;

//         for (auto& particle : particlesData(pk)) {
//             const double3 point_start = particle.coord;
//             const double3 velocity = particle.velocity;
//             double3 point_end = point_start + velocity * dt;
//             const double alpha = 0.5 * dt * charge / mass();
//             double3 velocity05;
//             bool covergence_cycle = false;
//             int iter = 0;
//             while (!covergence_cycle) {
//                 const double3 x05 = 0.5 * (point_start + point_end);
//                 std::cout << "interation " << iter << " " << point_start << " "
//                           << point_end << "\n";

//                 double length = (point_end - point_start).norm();
//                 double3 E = double3(0, 0, 0);
//                 double3 B = double3(0, 0, 0);
//                 bool is_move = true;
//                 double3 ray_start = point_start;
//                 double3 ray_end = point_end;
//                 while (is_move) {
//                     std::vector<int3> voxels =
//                         voxel_traversal(ray_start, ray_end, bin_size);
//                     double3 current_point = ray_start;
//                     double3 next_point;
//                     int3 final_voxel = voxels[voxels.size() - 1];

//                     make_periodic_bound(final_voxel);

//                     for (int i = 0; i < voxels.size(); i++) {
//                         if (i != voxels.size() - 1) {
//                             double t = find_ray_voxel_intersection_parameter(
//                                 ray_start, ray_end, voxels[i], voxels[i + 1],
//                                 bin_size);
//                             next_point =
//                                 get_point_in_ray(ray_start, ray_end, t);
//                         } else {
//                             next_point = ray_end;
//                         }
//                         Node node(current_point, bin_size);
//                         Node new_node(next_point, bin_size);
//                         fill_shape(node, shape, false);
//                         fill_shape(new_node, new_shape, false);
//                         double current_length =
//                             (current_point - next_point).norm();
//                         E += (current_length / length) *
//                              interpolateE_Chen(fieldE, node, shape, new_shape);
//                         B += (current_length / length) *
//                              interpolateB(fieldB, node, shape, new_shape);
//                         std::cout << "track: " << current_point << " "
//                                   << next_point << std::endl;
//                         // std::cout << E << std::endl;
//                         current_point = next_point;
//                         if (make_periodic_bound(voxels[i + 1], current_point,
//                                                 ray_end)) {
//                             ray_start = current_point;
//                             break;
//                         }
//                         is_move = (voxels[i] != final_voxel);
//                     }
//                 }

//                 double3 ap = velocity + alpha * E;
//                 velocity05 = (ap + alpha * cross(ap, B) +
//                               alpha * alpha * dot(ap, B) * B) /
//                              (1. + pow(alpha * B.norm(), 2));
//                 std::cout << "velocity05 " << velocity05 << std::endl;
//                 double3 xn = point_start + dt * velocity05;
//                 make_periodic_bound(ray_end);
//                 iter++;
//                 covergence_cycle = (xn - ray_end).norm() < 1e-8;
//                 ray_end = xn;
//             }
//             particle.coord = point_start + dt * velocity05;
//             particle.velocity = 2. * velocity05 - velocity;
//         }
//     }
// }

// bool ParticlesArray::track_particle(double3& coord, double3& velocity, const Field3d& fieldE,
//                                     const Field3d& fieldB, double dt) {
//                                         bool status = true;
//     double3 cellSize(xCellSize, yCellSize, zCellSize);
//     double bin_size = xCellSize;

//     ShapeK shape, new_shape;
//     const double3 ray_start = coord;
//     double3 ray_end = ray_start + velocity * dt;
//     const double alpha = 0.5 * dt * charge / mass();
//     double3 velocity05;
//     bool covergence_cycle = false;
//     int iter = 0;
//     // may be out if xn is not in area and need change dt
//     while (!covergence_cycle) {
//         const double3 x05 = 0.5 * (ray_start + ray_end);
//         std::cout << "interation " << iter << " " << ray_start << " " << ray_end
//                   << "\n";
//         std::vector<int3> voxels =
//             voxel_traversal(ray_start, ray_end, xCellSize);
//         double3 current_point = ray_start;
//         double3 next_point;
//         double length = (ray_end - ray_start).norm();
//         double3 E = double3(0, 0, 0);
//         double3 B = double3(0, 0, 0);
//         for (int i = 0; i < voxels.size(); i++) {
//             if (i != voxels.size() - 1) {
//                 double t = find_ray_voxel_intersection_parameter(
//                     ray_start, ray_end, voxels[i], voxels[i + 1], bin_size);
//                 next_point = get_point_in_ray(ray_start, ray_end, t);
//             } else {
//                 next_point = ray_end;
//             }
//             Node node(current_point, cellSize);
//             Node new_node(next_point, cellSize);
//             fill_shape(node, shape, false);
//             fill_shape(new_node, new_shape, false);
//             double current_length = (current_point - next_point).norm();
//             E += (current_length / length) *
//                  interpolateE_Chen(fieldE, node, shape, new_shape);
//             B += (current_length / length) *
//                  interpolateB(fieldB, node, shape, new_shape);
//             std::cout << "track: " << current_point << " " << next_point
//                       << std::endl;
//             // std::cout << E << std::endl;
//             current_point = next_point;
//         }
//         double3 ap = velocity + alpha * E;
//         velocity05 =
//             (ap + alpha * cross(ap, B) + alpha * alpha * dot(ap, B) * B) /
//             (1. + pow(alpha * B.norm(), 2));
//         std::cout << "velocity05 " << velocity05 << std::endl;
//         double3 xn = ray_start + dt * velocity05;
//         status = in_extended_domain(xn);
//         if (!status) {
//             return status;
//         }
//         iter++;
//         covergence_cycle = (xn - ray_end).norm() < 1e-8;
//         ray_end = xn;
//     }
//     coord = ray_end;
//     velocity = 2. * velocity05 - velocity;
//     return status;
// }

bool ParticlesArray::is_voxel_in_area(const int3& voxel) {
    int3 numCells(xCellCount, yCellCount, zCellCount);

    for (int dim = 0; dim < 3; ++dim) {
        if (voxel(dim) < 0 || voxel(dim) >= numCells(dim)) {
            return false;
        }
    }
    return true;
}

bool ParticlesArray::make_periodic_bound_force(double3& point) {
    bool make_bound = false;
    double eps = 1.e-13;
    double3 size = double3(static_cast<double>(xCellCount),
                           static_cast<double>(yCellCount),
                           static_cast<double>(zCellCount))*xCellSize;
    double3 new_point = point;
    for (int dim = 0; dim < 3; ++dim) {
        //std::cout << "point " << point(dim) << " " << size(dim) << std::endl;
        if (point(dim) <= 0) {
            new_point(dim) = point(dim) + size(dim) - eps;
            make_bound = true;
        }
        if ( size(dim) <= point(dim)) {
            new_point(dim) = point(dim) - size(dim) + eps;
            make_bound = true;
        }
    }
    point = new_point;
    //std::cout << "new point " << point << std::endl;
    return make_bound;
}

double ParticlesArray::track_particle(double3& coord, double3& velocity,
                                    const Field3d& fieldE,
                                    const Field3d& fieldB, double dt, bool& intersect_bound) {
    double3 cellSize(xCellSize, yCellSize, zCellSize);
    double bin_size = xCellSize;

    ShapeK shape, new_shape;
    const double3 ray_start = coord;
    double3 ray_end = ray_start + velocity * dt;
    double3 velocity05;
    double3 ray_end_inner; 
    bool covergence_cycle = false;
    int iter = 0;
    double dt1;
    // may be out if xn is not in area and need change dt
    while (!covergence_cycle) {
        dt1 = dt;
        const double3 x05 = 0.5 * (ray_start + ray_end);
        // std::cout << "interation " << iter << " " << ray_start << " " << ray_end
        //           << "\n";
        std::vector<int3> voxels =
            voxel_traversal(ray_start, ray_end, xCellSize);
        double3 current_point = ray_start;
        double3 next_point;
        ray_end_inner = ray_end;
        double3 E = double3(0, 0, 0);
        double3 B = double3(0, 0, 0);
        for (int i = 0; i < voxels.size(); i++) {
            if (!is_voxel_in_area(voxels[i])) {
                dt1 = dt * (ray_start - current_point).norm() /
                      (ray_start - ray_end).norm();
                ray_end_inner = current_point;
              //  std::cout << dt1 << " " << voxels[i] << " " << ray_start << " " << ray_end << "\n";
                break;
            }

            if (i != voxels.size() - 1) {
                double t = find_ray_voxel_intersection_parameter(
                    ray_start, ray_end, voxels[i], voxels[i + 1], bin_size);
                next_point = get_point_in_ray(ray_start, ray_end, t);
            } else {
                next_point = ray_end;
            }

            Node node(current_point, cellSize);
            Node new_node(next_point, cellSize);
            fill_shape(node, shape, false);
            fill_shape(new_node, new_shape, false);
            double current_length = (current_point - next_point).norm();
            E += current_length *
                 interpolateE_Chen(fieldE, node, shape, new_shape);
            // B.z() += current_length * (0.2 * (1 - 0.8 * (0.5*(current_point.y()+next_point.y())- 3.4)));
            // B.x() = B.y() = 0.;
            double3 mid_point = 0.5 * (current_point + next_point);
            node.set(mid_point, cellSize);
            fill_shape(node, shape, false);
            fill_shape(node, new_shape, true);
            B += current_length * interpolateB(fieldB, node, shape, new_shape);
            // std::cout << interpolateB(fieldB, node, shape, new_shape).z() -
            //           (0.2 * (1 - 0.8 * (mid_point.y() - 3.4))) << "\n";
            // std::cout << "track: " << current_point << " " << next_point
            //           << std::endl;
            //  std::cout << E << std::endl;

            current_point = next_point;
        }
        double length = (ray_end_inner - ray_start).norm();
        E/=length;
        B/=length;

        const double alpha = 0.5 * dt1 * charge / mass();

        double3 ap = velocity + alpha * E;
        velocity05 =
            (ap + alpha * cross(ap, B) + alpha * alpha * dot(ap, B) * B) /
            (1. + pow(alpha * B.norm(), 2));
        //std::cout << "velocity05 " << velocity05 << std::endl;
        double3 xn = ray_start + dt1 * velocity05;
        //double3 ray_end_inner = ray_start + dt1 * velocity05;

        iter++;
        // std::cout << "iter " << iter << " " << ray_start << " " << ray_end << " " << xn
        //           << " " << ray_end_inner << std::endl;
        covergence_cycle = (xn - ray_end_inner).norm() < 1e-10;
        //ray_end = xn;
        ray_end = ray_start + dt * velocity05;
    }

   std::cout << "interations " << iter << "\n";

    velocity = 2. * velocity05 - velocity;

    intersect_bound = dt1 < dt;

    coord = (intersect_bound) ? ray_end_inner : ray_end;

    return dt1;
}

void ParticlesArray::push_Chen(const Field3d& fieldE, const Field3d& fieldB,
                               double dt) {
    // const double dt_eps = 1.e-6;
    //  double3 cellSize(xCellSize, yCellSize, zCellSize);
    //  double bin_size = xCellSize;
#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        // ShapeK shape, new_shape;

        for (auto& particle : particlesData(pk)) {
            //std::cout << particle << "\n";
            double dtt = 0;
            double dt2 = dt;
            while (dtt < dt) {
                double3 coord = particle.coord;
                double3 velocity = particle.velocity;
                bool intersect_bound = false;
                double dt_fact = track_particle(coord, velocity, fieldE, fieldB,
                                                dt2, intersect_bound);
                std::cout << "track done " << dt_fact << "\n";
                if (intersect_bound) {
                    std::cout << "intersect bound" << particle <<  std::endl;
                    if (!make_periodic_bound_force(coord))
                        std::cout << "error in bound correction"
                                  << std::endl;
                }
                dtt += dt_fact;

                dt2 = dt - dtt;
               // if (intersect_bound) {std::cout<< dtt << " " << dt2 << "\n"; exit(0);}

                particle.coord = coord;
                particle.velocity = velocity;
                std::cout << particle << "\n";
            }
        }
    }
}

// void ParticlesArray::push_Chen(const Field3d& fieldE, const Field3d& fieldB,
//                                double dt) {
//     double3 cellSize(xCellSize, yCellSize, zCellSize);
//     double bin_size = xCellSize;
// #pragma omp parallel for
//     for (auto pk = 0; pk < size(); ++pk) {
//         ShapeK shape, new_shape;

//         for (auto& particle : particlesData(pk)) {
//             const double3 ray_start = particle.coord;
//             const double3 velocity = particle.velocity;
//             double3 ray_end = ray_start + velocity * dt;
//             const double alpha = 0.5 * dt * charge / mass();
//             double3 velocity05;
//             bool covergence_cycle = false;
//             int iter = 0;
//             while (!covergence_cycle) {
//                 const double3 x05 = 0.5 * (ray_start + ray_end);
//                 std::cout << "interation " << iter << " " << ray_start << " "
//                           << ray_end << "\n";
//                 std::vector<int3> voxels =
//                     voxel_traversal(ray_start, ray_end, xCellSize);
//                 double3 current_point = ray_start;
//                 double3 next_point;
//                 double length = (ray_end - ray_start).norm();
//                 double3 E = double3(0, 0, 0);
//                 double3 B = double3(0, 0, 0);
//                 for (int i = 0; i < voxels.size(); i++) {
//                     if (i != voxels.size() - 1) {
//                         double t = find_ray_voxel_intersection_parameter(
//                             ray_start, ray_end, voxels[i], voxels[i + 1],
//                             bin_size);
//                         next_point = get_point_in_ray(ray_start, ray_end, t);
//                     } else {
//                         next_point = ray_end;
//                     }
//                     Node node(current_point, cellSize);
//                     Node new_node(next_point, cellSize);
//                     fill_shape(node, shape, false);
//                     fill_shape(new_node, new_shape, false);
//                     double current_length = (current_point - next_point).norm();
//                     E += (current_length / length) *
//                          interpolateE_Chen(fieldE, node, shape, new_shape);
//                     B += (current_length / length) *
//                          interpolateB(fieldB, node, shape, new_shape);
//                     std::cout << "track: " << current_point << " " << next_point
//                               << std::endl;
//                     // std::cout << E << std::endl;
//                     current_point = next_point;
//                 }
//                 double3 ap = velocity + alpha * E;
//                 velocity05 = (ap + alpha * cross(ap, B) +
//                               alpha * alpha * dot(ap, B) * B) /
//                              (1. + pow(alpha * B.norm(), 2));
//                 std::cout << "velocity05 " << velocity05 << std::endl;
//                 double3 xn = ray_start + dt * velocity05;
//                 iter++;
//                 covergence_cycle = (xn - ray_end).norm() < 1e-8;
//                 ray_end = xn;
//             }
//             particle.coord = ray_end;
//             particle.velocity = 2. * velocity05 - velocity;
//         }
//     }
// }
