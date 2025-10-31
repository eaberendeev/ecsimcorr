/**
 * @file cross_section.h
 * @author Morozov O. P.
 * @brief Defines constants and enums for collision cross-section calculations.
 */
#ifndef CROSS_SECTION_HPP
#define CROSS_SECTION_HPP

#define A 1.53e-25
#define B 0.6
#define C 0.56

#define MC2_EV 511000.0
#define EV_TO_MC2 (1.0 / MC2_EV)
#define HYDROGEN_IONIZATION_EV 13.5984346

#define P (HYDROGEN_IONIZATION_EV * EV_TO_MC2)

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

#define E_ion P // Порог ионизации (mc^2)
#define B_param (10.0 * EV_TO_MC2) // Параметр модели Опала (~10 эВ)

enum class CollisionType { IONIZATION, CHARGE_EXCHANGE, NULL_COLLISION };

enum class ElectronIonizationModel { Approximation, Opal };

void set_electron_ionization_model(ElectronIonizationModel model);

ElectronIonizationModel get_electron_ionization_model();

#endif // CROSS_SECTION_HPP
