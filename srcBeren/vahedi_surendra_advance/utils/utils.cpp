// utils.cpp

#include "utils.hpp"
#include "Vec.h"
#include <cmath>

// Функция вычисления энергии частицы (E = 0.5 * m * v^2, нормированная на mc^2)
double compute_energy(const double3& velocity, double mass) {
    return 0.5 * mass * velocity.square();
}

double compute_velocity(double E, double mass) {
    return sqrt(2.*E/mass);
}

// Функция перевода скорости в систему центра масс
double3 v_center_of_mass(double3 v1, double3 v2, double m1, double m2) {
    return (m1 * v1 + m2 * v2) / (m1 + m2);
}

// Функция перевода скорости обратно в лабораторную систему
double3 to_lab_frame(double3 v_com, double3 v_other, double m_self, double m_other) {
    return v_com - (m_other / (m_self + m_other)) * v_other;
}
