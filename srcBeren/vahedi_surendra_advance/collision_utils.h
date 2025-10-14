/**
 * @file collision_utils.h
 * @author Morozov O. P.
 * @brief Declares utility functions used when computing collision dynamics.
 */
#pragma once
#include "Vec.h"

double compute_energy(const double3& velocity, double mass);

double compute_velocity(double E, double mass);

double3 v_center_of_mass(double3 v1, double3 v2, double m1, double m2);

double3 to_lab_frame(double3 v_com, double3 v_other, double m_self, double m_other);

inline double compute_collision_probability(double freq, double dt) {
    return 1.0 - exp(-freq * dt);
}
