// utils.hpp

#ifndef UTILS_HPP
#define UTILS_HPP

//#include "clsVec.hpp"
#include "cross_section.hpp"
#include "Vec.h"

// Функция для вычисления энергии частицы
double compute_energy(const double3& velocity, double mass);

double compute_velocity(double E, double mass);

// Функция перевода скорости в систему центра масс
double3 v_center_of_mass(double3 v1, double3 v2, double m1, double m2);

// Функция перевода скорости обратно в лабораторную систему
double3 to_lab_frame(double3 v_com, double3 v_other, double m_self, double m_other);

#endif // UTILS_HPP
