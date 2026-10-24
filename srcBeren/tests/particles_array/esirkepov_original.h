#pragma once

#include "Particle.h"
#include "ParticlesArray.h"
#include "World.h"
#include "config.h"
#include "particle_path.h"
#include "timer.h"

// For understand algorithm. Slow version
template <ShapeFunction ShapeFn, int ShapeSize>
inline void move_and_calc_current_impl_original(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                                double charge, double mpw, const double dt, Field3d& fieldJ,
                                                const PathPolicy& policy = PathPolicy::none()) {
    RECORD_TIMER;
    constexpr auto SMAX = 2 * ShapeSize;

    const double conx = domain.cell_size().x() / (6 * dt) * mpw;
    const double cony = domain.cell_size().y() / (6 * dt) * mpw;
    const double conz = domain.cell_size().z() / (6 * dt) * mpw;

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < particlesData.capacity(); ++pk) {
        double arg;
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        for (auto& particle : particlesData(pk)) {
            Vector3R start = particle.coord;

            particle.move(dt);

            Vector3R end = particle.coord;
            Vector3R end_for_J = end;
            if (policy.active())
                shape_particle_path(start, end, domain, policy, end_for_J);

            double xx = start.x() / domain.cell_size().x();
            double yy = start.y() / domain.cell_size().y();
            double zz = start.z() / domain.cell_size().z();

            double xn = end_for_J.x() / domain.cell_size().x();
            double yn = end_for_J.y() / domain.cell_size().y();
            double zn = end_for_J.z() / domain.cell_size().z();

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
                sx[n] = ShapeFn(arg);
                arg = -yy + double(yk - GHOST_CELLS + n);
                sy[n] = ShapeFn(arg);
                arg = -zz + double(zk - GHOST_CELLS + n);
                sz[n] = ShapeFn(arg);
                arg = -xn + double(xk - GHOST_CELLS + n);
                sx_n[n] = ShapeFn(arg);
                arg = -yn + double(yk - GHOST_CELLS + n);
                sy_n[n] = ShapeFn(arg);
                arg = -zn + double(zk - GHOST_CELLS + n);
                sz_n[n] = ShapeFn(arg);
            }

            for (int n = 0; n < SMAX; ++n) {
                int indx = xk + n;
                for (int m = 0; m < SMAX; ++m) {
                    int indy = yk + m;
                    for (int k = 0; k < SMAX; ++k) {
                        int indz = zk + k;

                        if (n == 0)
                            jx[n][m][k] = -charge * conx * (sx_n[n] - sx[n]) *
                                          (sy_n[m] * (2 * sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                        if (n > 0 && n < SMAX - 1)
                            jx[n][m][k] =
                                jx[n - 1][m][k] - charge * conx * (sx_n[n] - sx[n]) *
                                                      (sy_n[m] * (2 * sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                        if (m == 0)
                            jy[n][m][k] = -charge * cony * (sy_n[m] - sy[m]) *
                                          (sx_n[n] * (2 * sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                        if (m > 0 && m < SMAX - 1)
                            jy[n][m][k] =
                                jy[n][m - 1][k] - charge * cony * (sy_n[m] - sy[m]) *
                                                      (sx_n[n] * (2 * sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));

                        if (k == 0)
                            jz[n][m][k] = -charge * conz * (sz_n[k] - sz[k]) *
                                          (sy_n[m] * (2 * sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                        if (k > 0 && k < SMAX - 1)
                            jz[n][m][k] =
                                jz[n][m][k - 1] - charge * conz * (sz_n[k] - sz[k]) *
                                                      (sy_n[m] * (2 * sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));

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

template <ShapeType ST = SHAPE>
inline void move_and_calc_current_original(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                           double charge, double mpw, const double dt, Field3d& fieldJ,
                                           const PathPolicy& policy = PathPolicy::none()) {
    RECORD_TIMER;
    static_assert(ST == ShapeType::Linear || ST == ShapeType::Quadratic,
                  "shape dispatch supports Linear and Quadratic only");
    if constexpr (ST == ShapeType::Linear) {
        move_and_calc_current_impl_original<Shape, 2>(particlesData, domain, charge, mpw, dt, fieldJ, policy);
    } else {
        move_and_calc_current_impl_original<Shape2, 2>(particlesData, domain, charge, mpw, dt, fieldJ, policy);
    }
}

template <ShapeType ST = SHAPE>
inline void move_and_calc_current_original(ParticlesArray& sp, double dt, Field3d& fieldJ,
                                           const PathPolicy& policy = PathPolicy::none()) {
    if (sp.is_neutral())
        return;
    ::move_and_calc_current_original<ST>(sp.particlesData, sp.get_domain(), sp.charge, sp.mpw(), dt, fieldJ, policy);
}
