#ifndef CROSS_SECTION_HPP
#define CROSS_SECTION_HPP

constexpr double CROSS_SECTION_A = 1.53e-25;
constexpr double CROSS_SECTION_B = 0.6;
constexpr double CROSS_SECTION_C = 0.56;

constexpr double MC2_EV = 511000.0;
constexpr double EV_TO_MC2 = 1.0 / MC2_EV;
constexpr double HYDROGEN_IONIZATION_EV = 13.5984346;

constexpr double THRESHOLD_IONIZATION = HYDROGEN_IONIZATION_EV * EV_TO_MC2;

constexpr double CROSS_SECTION_amuH = 1.00783;

constexpr double CROSS_SECTION_A1 = 3.2345;
constexpr double CROSS_SECTION_A2 = 235.88;
constexpr double CROSS_SECTION_A3 = 0.038371;
constexpr double CROSS_SECTION_A4 = 3.8068e-6;
constexpr double CROSS_SECTION_A5 = 1.1832e-10;
constexpr double CROSS_SECTION_A6 = 2.3713;

constexpr double CROSS_SECTION_B1 = 12.899;
constexpr double CROSS_SECTION_B2 = 61.897;
constexpr double CROSS_SECTION_B3 = 9.2731e3;
constexpr double CROSS_SECTION_B4 = 4.9749e-4;
constexpr double CROSS_SECTION_B5 = 3.9890e-2;
constexpr double CROSS_SECTION_B6 = -1.5900;
constexpr double CROSS_SECTION_B7 = 3.1834;
constexpr double CROSS_SECTION_B8 = -3.7154;

#define E_ion   THRESHOLD_IONIZATION   // Порог ионизации (mc^2)
#define B_param (10.0 * EV_TO_MC2)     // Параметр модели Опала (~10 эВ)

enum class CollisionType { IONIZATION, CHARGE_EXCHANGE, NULL_COLLISION };

enum class ElectronIonizationModel { Approximation, Opal };

void set_electron_ionization_model(ElectronIonizationModel model);

ElectronIonizationModel get_electron_ionization_model();

#endif   // CROSS_SECTION_HPP
