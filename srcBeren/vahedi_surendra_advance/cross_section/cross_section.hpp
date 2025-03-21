// cross_section.hpp

#ifndef CROSS_SECTION_HPP
#define CROSS_SECTION_HPP

#define mc2 511. // keV
#define c 3e10 // см/c
#define Me 1. // в массах электрона
#define Mi 1836. // в массах электрона
#define Mn 1836. // в массах электрона

#define A 3.82964e-26
#define B 0.6
#define C 0.56
#define P 2.6614481409001956e-05 // безразмерная энергия

#define amuH 1.00783

#define A1 3.2345
#define A2 235.88
#define A3 0.038371
#define A4 3.8068e-6
#define A5 1.1832e-10
#define A6 2.3713

#define B1 12.899
#define B2 61.897
#define B3 9.2731e3
#define B4 4.9749e-4
#define B5 3.9890e-2
#define B6 -1.5900
#define B7 3.1834
#define B8 -3.7154

#define E_ion 0. // Порог ионизации
#define B_param 10. // Параметр модели Опала (порядка 10 эВ)

#define l0 0.001 // с/w_pl
#define dt 0.001
#define n0 1e13

// Функции для вычисления сечений столкновений
double Sigma_e(double E); // Ударная ионизация электронами
double Sigma_p(double E); // Ударная ионизация протонами
double Sigma_cx(double E); // Резонансная перезарядка

#endif // CROSS_SECTION_HPP
