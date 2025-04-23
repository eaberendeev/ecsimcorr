#include <functional>

#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "service.h"
#include "sgs.h"

using ShapeFunction = double (*)(const double&);

// Template version for compile-time shape function selection
template <ShapeFunction ShapeFn, int ShapeSize>
void ParticlesArray::density_on_grid_update_impl() {
    constexpr auto SMAX = 2 * ShapeSize;
    densityOnGrid.setZero();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto j = 0; j < size(); ++j) {
        alignas(64) double sx[SMAX];
        alignas(64) double sy[SMAX];
        alignas(64) double sz[SMAX];

        for (const auto& particle : particlesData(j)) {
            // Vectorizable coordinate calculations
            const double ix = particle.coord.x() / xCellSize;
            const double iy = particle.coord.y() / yCellSize;
            const double iz = particle.coord.z() / zCellSize;

            const int xk = int(ix);
            const int yk = int(iy);
            const int zk = int(iz);

// Vectorizable shape calculations
#pragma omp simd
            for (int n = 0; n < SMAX; ++n) {
                sx[n] = ShapeFn(-ix + double(xk - GHOST_CELLS + n));
                sy[n] = ShapeFn(-iy + double(yk - GHOST_CELLS + n));
                sz[n] = ShapeFn(-iz + double(zk - GHOST_CELLS + n));
            }

            const double weight = _mpw * charge;

// Density accumulation with loop unrolling
#pragma unroll
            for (int n = 0; n < SMAX; ++n) {
                const int indx = xk + n;
                const double sxw = sx[n] * weight;

                for (int m = 0; m < SMAX; ++m) {
                    const int indy = yk + m;
                    const double sxyw = sxw * sy[m];

                    for (int k = 0; k < SMAX; ++k) {
                        const int indz = zk + k;
#pragma omp atomic update
                        densityOnGrid(indx, indy, indz, 0) += sxyw * sz[k];
                    }
                }
            }
        }
    }
    apply_periodic_border_with_add(densityOnGrid, bounds);
}

template <typename VelocityCalculator1, typename VelocityCalculator2>
void ParticlesArray::calculate_pressure_component(
    Field3d& P, VelocityCalculator1 velocityCalc1,
    VelocityCalculator2 velocityCalc2) {
    P.setZero();
    constexpr auto SMAX = 2;
    const double x0 = 0.5 * xCellSize * xCellCount;

#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < size(); ++pk) {
        double wx[SMAX], wy[SMAX], wz[SMAX];

        for (auto& particle : particlesData(pk)) {
            const auto& coord = particle.coord;
            const auto& velocity = particle.velocity;

            auto x = coord.x() / xCellSize + GHOST_CELLS;
            auto y = coord.y() / yCellSize + GHOST_CELLS;
            auto z = coord.z() / zCellSize + GHOST_CELLS;

            const auto intx = int(x);
            const auto inty = int(y);
            const auto intz = int(z);

            wx[1] = (x - intx);
            wx[0] = 1 - wx[1];
            wy[1] = (y - inty);
            wy[0] = 1 - wy[1];
            wz[1] = (z - intz);
            wz[0] = 1 - wz[1];

            double v1 = velocityCalc1(coord, velocity, x0);
            double v2 = velocityCalc2(coord, velocity, x0);
            double pressure = _mass * v1 * v2 * _mpw;

            for (int nx = 0; nx < SMAX; ++nx) {
                const int i = intx + nx;
                for (int ny = 0; ny < SMAX; ++ny) {
                    const int j = inty + ny;
                    for (int nz = 0; nz < SMAX; ++nz) {
                        const int k = intz + nz;
                        const auto sx = wx[nx] * wy[ny] * wz[nz];
#pragma omp atomic update
                        P(i, j, k, 0) += sx * pressure;
                    }
                }
            }
        }
    }
}
