#include "cross_section.h"
#include "collisions_with_neutrals.h"
#include <cmath>

// Функции сечений столкновений

namespace {
ElectronIonizationModel g_electron_ionization_model = ElectronIonizationModel::Approximation;
constexpr double kIonizationConsistency = (THRESHOLD_IONIZATION > E_ion) ? (THRESHOLD_IONIZATION - E_ion) : (E_ion - THRESHOLD_IONIZATION);
static_assert(kIonizationConsistency < 1e-12, "E_ion должно совпадать с THRESHOLD_IONIZATION в выбранных единицах");
}

void set_electron_ionization_model(ElectronIonizationModel model) {
    g_electron_ionization_model = model;
}

ElectronIonizationModel get_electron_ionization_model() {
    return g_electron_ionization_model;
}

// Ионизация электронным ударом:
double ColliderWithNeutrals::Sigma_e(double energy) const {
    if (energy > THRESHOLD_IONIZATION) {
        double s1 = CROSS_SECTION_B *
                    exp(-CROSS_SECTION_C * (energy / THRESHOLD_IONIZATION - 1.));
        double s2 = CROSS_SECTION_A * log(energy / THRESHOLD_IONIZATION) /
                    (energy * THRESHOLD_IONIZATION);
        return s2*(1 - s1)*dpd0;
    } else {return 0.;}
}

// Резонансная перезарядка:
double ColliderWithNeutrals::Sigma_cx(double energy) const {
    energy = energy * SGS::MC2 / CROSS_SECTION_amuH;
    double Scx = 1e-16;
    Scx *= (CROSS_SECTION_A1 * log(CROSS_SECTION_A2 / energy + CROSS_SECTION_A6));
    Scx /= (1 + CROSS_SECTION_A3 * energy + CROSS_SECTION_A4 * pow(energy, 3.5) +
            CROSS_SECTION_A5 * pow(energy, 5.4));
    return Scx*dpd0;
}

// Резонансная перезарядка:
double ColliderWithNeutrals::Sigma_p(double energy) const {
    if (energy > THRESHOLD_IONIZATION) {
        energy = energy * SGS::MC2 / CROSS_SECTION_amuH;
        double Sp = 1e-16;
        double Sp_left =
            exp(-CROSS_SECTION_B2 / energy) * log(1. + CROSS_SECTION_B3 * energy) / energy;
        double Sp_right = CROSS_SECTION_B4 * exp(-CROSS_SECTION_B5 * energy) /
                          (pow(energy, CROSS_SECTION_B6) +
                           CROSS_SECTION_B7 * pow(energy, CROSS_SECTION_B8));
        Sp *= CROSS_SECTION_B1 * (Sp_left + Sp_right);
        return Sp*dpd0;
    } else {return 0.;}
}
