#pragma once

#include "Particle.h"
#include "vector3.h"

namespace borisPusher {

inline void update_vEB(Particle& particle, const double qm, const Vector3R& E_p, const Vector3R& B_p, const double dt) {
    const double alpha = dt * qm;
    const Vector3R a = +alpha * E_p;
    const Vector3R b = -alpha * B_p;

    Vector3R& v = particle.velocity;
    const Vector3R w = v + 0.5 * a;
    const Vector3R tmp = b.cross(w);
    v += a + (tmp + 0.5 * b.cross(tmp)) / (1.0 + 0.25 * b.squared());
}

}   // namespace borisPusher
