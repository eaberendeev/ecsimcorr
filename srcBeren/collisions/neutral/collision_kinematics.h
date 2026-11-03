#pragma once

#include <cmath>

#include "core/vector3.h"

// Energy and speed use the particle code normalization, without conversion to eV.
inline double compute_energy(const Vector3R& velocity, double mass) {
    return 0.5 * mass * velocity.squaredNorm();
}

inline double compute_velocity(double E, double mass) {
    return sqrt(2. * E / mass);
}

inline Vector3R v_center_of_mass(Vector3R v1, Vector3R v2, double m1, double m2) {
    return (m1 * v1 + m2 * v2) / (m1 + m2);
}

inline double compute_collision_probability(double freq, double dt) {
    return 1.0 - exp(-freq * dt);
}
