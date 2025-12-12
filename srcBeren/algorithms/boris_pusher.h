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
    v += a + (b.cross(w) + 0.5 * b.cross(b.cross(w))) /
                 (1.0 + 0.25 * b.squared());
}

}   // namespace BorisPusher