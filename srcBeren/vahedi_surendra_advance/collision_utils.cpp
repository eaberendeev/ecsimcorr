#include "collision_utils.h"

#include <cmath>

double compute_energy(const Vector3R& velocity, double mass) {
    return 0.5 * mass * velocity.squared();
}

double compute_velocity(double E, double mass) {
    return sqrt(2. * E / mass);
}

Vector3R v_center_of_mass(Vector3R v1, Vector3R v2, double m1, double m2) {
    return (m1 * v1 + m2 * v2) / (m1 + m2);
}

Vector3R to_lab_frame(Vector3R v_com, Vector3R v_other, double m_self, double m_other) {
    return v_com - (m_other / (m_self + m_other)) * v_other;
}
