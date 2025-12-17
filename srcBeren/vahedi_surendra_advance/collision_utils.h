#pragma once
#include "vector3.h"

double compute_energy(const Vector3R& velocity, double mass);

double compute_velocity(double E, double mass);

Vector3R v_center_of_mass(Vector3R v1, Vector3R v2, double m1, double m2);

Vector3R to_lab_frame(Vector3R v_com, Vector3R v_other, double m_self, double m_other);

inline double compute_collision_probability(double freq, double dt) {
    return 1.0 - exp(-freq * dt);
}
