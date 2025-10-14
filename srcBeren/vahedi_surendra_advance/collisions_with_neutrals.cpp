/**
 * @file collisions_with_neutrals.cpp
 * @author Morozov O. P.
 * @brief Implements neutral collision calculations and data export helpers.
 */
#include "collisions_with_neutrals.h"

#include "collision_processing.h"
#include "collision_utils.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

std::pair<double, double> ColliderWithNeutrals::compute_frequencies(const double3& vcp,
                                                                   const double3& vn,
                                                                   double mcp,
                                                                   double nn) const {
    bool is_electron = (mcp == 1.0);
    bool is_proton = (mcp == 1836.0);

    if (!is_electron && !is_proton) {
        std::cout << " mcp " << mcp << "\n";
        throw std::invalid_argument("Неизвестный тип заряженной частицы");
    }

    double3 v_rel = vcp - vn;
    double E = compute_energy(v_rel, mcp);
    double v_mod = compute_velocity(E, mcp);

    double ion_freq = 0.0;
    double cx_freq = 0.0;

    if (is_electron) {
        if (process_options.electron_ionization) {
            ion_freq = v_mod * nn * Sigma_e(E);
        }
    } else {
        if (process_options.proton_ionization) {
            ion_freq = v_mod * nn * Sigma_p(E);
        }
        if (process_options.proton_charge_exchange) {
            cx_freq = v_mod * nn * Sigma_cx(E);
        }
    }


    return {ion_freq, cx_freq};
}

double ColliderWithNeutrals::total_collision_frequency(const double3& vcp,
                                                       const double3& vn,
                                                       double mcp,
                                                       double nn) const {
    auto [ion_freq, cx_freq] = compute_frequencies(vcp, vn, mcp, nn);
    return ion_freq + cx_freq;
}

std::tuple<bool, double3, double3> ColliderWithNeutrals::collision_with_neutral(
    double3& vcp, double3& vn, double mcp, double mn, double ncp, double nn, double dt,
    double freq_max
) {
    bool is_electron = (mcp == 1.0);
    bool is_proton = (mcp == 1836.0);

    if (!is_electron && !is_proton) {
        std::cout << " mcp " << mcp << "\n";
        throw std::invalid_argument("Неизвестный тип заряженной частицы");
    }

    if (scheme_mode == CollisionScheme::NULL_COLLISION) {
        double P_null_collision = compute_collision_probability(freq_max, dt);
        if (!check_collision(P_null_collision)) {
            return {false, vcp, vn};
        }
        auto [ion_freq, cx_freq] = compute_frequencies(vcp, vn, mcp, nn);
        double effective_freq_max = std::max(freq_max, ion_freq + cx_freq);
        CollisionType collision_type = select_collision_type(is_electron, ion_freq, cx_freq, effective_freq_max);
        return process_collision(collision_type, is_electron, vcp, vn, mcp, mn);
    }

    auto [ion_freq, cx_freq] = compute_frequencies(vcp, vn, mcp, nn);

    double total_freq = ion_freq + cx_freq;
    if (total_freq <= 0.0) {
        return {false, vcp, vn};
    }

    double P_collision = compute_collision_probability(total_freq, dt);
    if (!check_collision(P_collision)) {
        return {false, vcp, vn};
    }

    CollisionType collision_type = select_collision_type(is_electron, ion_freq, cx_freq, total_freq);
    return process_collision(collision_type, is_electron, vcp, vn, mcp, mn);
}
