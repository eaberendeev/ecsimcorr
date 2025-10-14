/**
 * @file cross_section.cpp
 * @author Morozov O. P.
 * @brief Implements cross-section models for particle-neutral interactions.
 */
#include "cross_section.h"
#include "collisions_with_neutrals.h"
#include <cmath>

// Функции сечений столкновений

namespace {
ElectronIonizationModel g_electron_ionization_model = ElectronIonizationModel::Approximation;
constexpr double kIonizationConsistency = (P > E_ion) ? (P - E_ion) : (E_ion - P);
static_assert(kIonizationConsistency < 1e-12, "E_ion должно совпадать с P в выбранных единицах");
}

void set_electron_ionization_model(ElectronIonizationModel model) {
    g_electron_ionization_model = model;
}

ElectronIonizationModel get_electron_ionization_model() {
    return g_electron_ionization_model;
}

// Ионизация электронным ударом:
double ColliderWithNeutrals::Sigma_e(double E) const {
    if (E > P) {
        double s1 = B * exp(-C * (E/P - 1.));
        double s2 = A * log(E/P) / (E*P);
        return s2*(1 - s1)*dpd0;
    } else {return 0.;}
}

// Резонансная перезарядка:
double ColliderWithNeutrals::Sigma_cx(double E) const {
    E = E*SGS::MC2/amuH;
    double Scx = 1e-16;
    Scx *= (A1 * log(A2/E + A6));
    Scx /= (1+A3*E+A4*pow(E, 3.5)+A5*pow(E, 5.4));
    return Scx*dpd0;
}

// Резонансная перезарядка:
double ColliderWithNeutrals::Sigma_p(double E) const {
    if (E > P) {
        E = E*SGS::MC2/amuH;
        double Sp = 1e-16;
        double Sp_left = exp(- B2 /E) * log(1.+ B3 *E) / E;
        double Sp_right = B4 * exp(-B5 * E) / (pow(E, B6) + B7 * pow(E, B8)); 
        Sp *= B1 * (Sp_left + Sp_right);
        return Sp*dpd0;
    } else {return 0.;}
}
