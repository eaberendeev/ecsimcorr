#pragma once

#include "Particle.h"
#include "Vec.h"

namespace borisPusher {

inline void update_vEB(Particle& particle, const double qm, const double3& E_p,
                       const double3& B_p, const double dt) {
    double alpha = dt * qm;
    double3 a = +alpha * E_p;
    double3 b = -alpha * B_p;

    double3& v = particle.velocity;
    double3 w = v + 0.5 * a;
    v += a + (cross(b, w) + 0.5 * cross(b, cross(b, w))) /
                 (1.0 + 0.25 * b.square());
}

}   // namespace BorisPusher