#pragma once

#include <functional>

#include "ParticlesArray.h"
#include "coupling/shape.h"
#include "collision.h"
#include "containers.h"
#include "grid/domain.h"
#include "sgs.h"

template <ShapeType ST>
inline void ParticlesArray::density_on_grid_update(Field3d& density) const {
    density_on_grid_update_impl<ShapeKernel<ST>::weight, ShapeKernel<ST>::deposit_half_width>(density);
}

template <ShapeType ST>
inline void ParticlesArray::density_on_grid_update_reference(Field3d& density) const {
    density_on_grid_update_impl_reference<ShapeKernel<ST>::weight, ShapeKernel<ST>::deposit_half_width>(density);
}

template <typename VelocityCalculator1, typename VelocityCalculator2>
void ParticlesArray::calculate_pressure_component(Field3d& P, VelocityCalculator1 velocityCalc1,
                                                  VelocityCalculator2 velocityCalc2) const {
    RECORD_TIMER;
    P.setZero();
    constexpr auto SMAX = 2;
    const Vector3R center =
        domain_.geom.use_cylinder ? domain_.geom.cyl_center : 0.5 * (domain_.geom.box_min + domain_.geom.box_max);

#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < size(); ++pk) {
        double wx[SMAX], wy[SMAX], wz[SMAX];

        for (auto& particle : particlesData(pk)) {
            const auto& coord = particle.coord;
            const auto& velocity = particle.velocity;

            auto x = coord.x() / domain_.cell_size().x() + GHOST_CELLS;
            auto y = coord.y() / domain_.cell_size().y() + GHOST_CELLS;
            auto z = coord.z() / domain_.cell_size().z() + GHOST_CELLS;

            const auto intx = int(x);
            const auto inty = int(y);
            const auto intz = int(z);

            wx[1] = (x - intx);
            wx[0] = 1 - wx[1];
            wy[1] = (y - inty);
            wy[0] = 1 - wy[1];
            wz[1] = (z - intz);
            wz[0] = 1 - wz[1];

            double v1 = velocityCalc1(coord, velocity, center);
            double v2 = velocityCalc2(coord, velocity, center);
            double pressure = mass_ * v1 * v2 * mpw_;

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
