#pragma once

#include "Particle.h"
#include "vector3.h"

namespace borisPusher {

inline void update_vEB(Particle& particle, const double qm, const Vector3R& E_p,
                       const Vector3R& B_p, const double dt) {
    double alpha = dt * qm;
    Vector3R a = +alpha * E_p;
    Vector3R b = -alpha * B_p;

    Vector3R& v = particle.velocity;
    Vector3R w = v + 0.5 * a;
    v += a + (b.cross(w) + 0.5 * b.cross(b.cross(w))) /
                 (1.0 + 0.25 * b.squared());
}

}   // namespace BorisPusher