/**
 * @file collisions_with_neutrals.h
 * @author Morozov O. P.
 * @brief Declares the neutral collision handler and related configuration types.
 */
#pragma once

#include <tuple>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <iostream>

#include "Vec.h"
#include "cross_section.h"
#include "sgs.h"

enum class CollisionScheme {
    NULL_COLLISION,
    PHYSICAL_ONLY
};

struct CollisionProcessOptions {
    bool electron_ionization = false;
    bool proton_ionization = false;
    bool proton_charge_exchange = true;
};

class ColliderWithNeutrals {
   public:
    explicit ColliderWithNeutrals(const double _n0,
                                  CollisionScheme scheme = CollisionScheme::PHYSICAL_ONLY,
                                  CollisionProcessOptions process_opts = CollisionProcessOptions()) {
        n0 = _n0;
        wp = SGS::get_plasma_freq(_n0);
        l0 = SGS::c / wp;
        dpd0 = n0 * l0;
        scheme_mode = scheme;
        process_options = process_opts;
        std::cout << "ColliderWithNeutrals create with parametrs: \n";
        std::cout << "n0: " << _n0 << "\n";
        std::cout << "wp: " << wp << "\n";
        std::cout << "l0: " << l0 << "\n"; 
        std::cout << "dpd0: " << dpd0 << "\n";
    }

    void set_scheme(CollisionScheme new_scheme) { scheme_mode = new_scheme; }
    CollisionScheme get_scheme() const { return scheme_mode; }
    void set_process_options(const CollisionProcessOptions& new_options) { process_options = new_options; }
    const CollisionProcessOptions& get_process_options() const { return process_options; }
    std::tuple<bool, double3, double3> collision_with_neutral(
        double3& vcp, double3& vn, double mcp, double mn, double ncp, double nn, double dt,
        double freq_max
    );

    double total_collision_frequency(const double3& vcp,
                                     const double3& vn,
                                     double mcp,
                                     double nn) const;

   private:
    double n0;
    double wp;
    double l0;
    double dpd0;

    double Sigma_e(double E) const;  
    double Sigma_p(double E) const;
    double Sigma_cx(double E) const; 

    CollisionProcessOptions process_options;

    std::pair<double, double> compute_frequencies(const double3& vcp,
                                                  const double3& vn,
                                                  double mcp,
                                                  double nn) const;

    bool check_collision(double P_collision);

    CollisionType select_collision_type(bool is_electron, double ion_freq, double cx_freq,
                                        double freq_bound);

    CollisionScheme scheme_mode;
};
