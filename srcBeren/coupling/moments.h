#pragma once

#include "ParticlesArray.h"
#include "timer.h"

namespace coupling {
// The cylindrical basis is undefined exactly on the axis; choose the x/y
// basis there so moments remain finite and the pressure trace is preserved.
inline Vector3R cylindrical_velocity(const Vector3R& coord, const Vector3R& velocity, const Vector3R& center) {
    const double x = coord.x() - center.x(), y = coord.y() - center.y();
    const double radius2 = x * x + y * y;
    if (radius2 == 0)
        return velocity;
    const double inverse_radius = 1.0 / std::sqrt(radius2);
    return {(x * velocity.x() + y * velocity.y()) * inverse_radius,
            (-y * velocity.x() + x * velocity.y()) * inverse_radius, velocity.z()};
}

inline double radial_velocity(const Vector3R& coord, const Vector3R& velocity, const Vector3R& center) {
    return cylindrical_velocity(coord, velocity, center).x();
}
inline double phi_velocity(const Vector3R& coord, const Vector3R& velocity, const Vector3R& center) {
    return cylindrical_velocity(coord, velocity, center).y();
}
inline double z_velocity(const Vector3R&, const Vector3R& velocity, const Vector3R&) {
    return velocity.z();
}

// Overwrite P with the raw second velocity moment m*mpw*sum(v1*v2*S),
// using the linear scalar-node shape. The mean flow is not subtracted:
// this is not thermal pressure. Coordinates and velocities are at the caller's
// supplied particle time level; the caller applies any required boundary BC.
template <typename VelocityCalculator1, typename VelocityCalculator2>
void deposit_second_moment(const ParticlesArray& species, Field3d& P, VelocityCalculator1 velocityCalc1,
                           VelocityCalculator2 velocityCalc2) {
    RECORD_TIMER;
    const Array3D<std::vector<Particle>>& particlesData = species.particlesData;
    const Domain& domain_ = species.get_domain();
    const double mass_ = species.mass();
    const double mpw_ = species.mpw();
    P.setZero();
    constexpr int SMAX = 2;
    const Vector3R center =
        domain_.geom.use_cylinder ? domain_.geom.cyl_center : 0.5 * (domain_.geom.box_min + domain_.geom.box_max);

#pragma omp parallel for schedule(dynamic, species.omp_granularity())
    for (int pk = 0; pk < species.size(); ++pk) {
        double wx[SMAX], wy[SMAX], wz[SMAX];

        for (const Particle& particle : particlesData(pk)) {
            const Vector3R& coord = particle.coord;
            const Vector3R& velocity = particle.velocity;

            double x = coord.x() / domain_.cell_size().x() + GHOST_CELLS;
            double y = coord.y() / domain_.cell_size().y() + GHOST_CELLS;
            double z = coord.z() / domain_.cell_size().z() + GHOST_CELLS;

            const int intx = int(x);
            const int inty = int(y);
            const int intz = int(z);

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
                        const double sx = wx[nx] * wy[ny] * wz[nz];
#pragma omp atomic update
                        P(i, j, k, 0) += sx * pressure;
                    }
                }
            }
        }
    }
}

}   // namespace coupling
