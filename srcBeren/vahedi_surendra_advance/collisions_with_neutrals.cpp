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
#include <omp.h>

std::pair<double, double> ColliderWithNeutrals::compute_frequencies(const double3& vcp,
                                                                   const double3& vn,
                                                                   double mcp,
                                                                   double nn) const {
    bool is_electron = (mcp  < 2);
    bool is_proton = (mcp >= 2);

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
    double3& vcp, double3& vn, double mcp, double mn, double ncp, double nn,
    double dt, double freq_max) {
    using clock = std::chrono::high_resolution_clock;
    auto t_start_total = clock::now();
    profiler.calls += 1;

    bool is_electron = (mcp < 2);
    bool is_proton = (mcp >= 2);

    if (!is_electron && !is_proton) {
        std::cout << " mcp " << mcp << "\n";
        throw std::invalid_argument("Неизвестный тип заряженной частицы");
    }

    // ----------------- NULL_COLLISION scheme -----------------
    if (scheme_mode == CollisionScheme::NULL_COLLISION) {
        // compute probability for null collision (timed)
        double P_null_collision;
        {
            auto t_cp_s = clock::now();
            P_null_collision = compute_collision_probability(freq_max, dt);
            auto t_cp_e = clock::now();
            profiler.time_compute_prob_ns +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(t_cp_e -
                                                                     t_cp_s)
                    .count();
        }

        // check collision (timed)
        bool happened;
        {
            auto t_check_s = clock::now();
            happened = check_collision(P_null_collision);
            auto t_check_e = clock::now();
            profiler.time_check_collision_ns +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(t_check_e -
                                                                     t_check_s)
                    .count();
        }

        if (!happened) {
            auto t_end_total = clock::now();
            profiler.time_total_ns +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t_end_total - t_start_total)
                    .count();
            return {false, vcp, vn};
        }

        // If happened -> compute frequencies (timed)
        double ion_freq = 0.0, cx_freq = 0.0;
        {
            auto t_cf0 = clock::now();
            std::tie(ion_freq, cx_freq) = compute_frequencies(vcp, vn, mcp, nn);
            auto t_cf1 = clock::now();
            profiler.time_compute_freq_ns +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(t_cf1 -
                                                                     t_cf0)
                    .count();
        }

        // record sample: P used here = P_null_collision (probability that was
        // checked)
        profiler.add_freq_sample(ion_freq, cx_freq, P_null_collision);

        // select collision type (timed)
        CollisionType collision_type;
        {
            auto t_sel_s = clock::now();
            double effective_freq_max = std::max(freq_max, ion_freq + cx_freq);
            collision_type = select_collision_type(is_electron, ion_freq,
                                                   cx_freq, effective_freq_max);
            auto t_sel_e = clock::now();
            profiler.time_select_type_ns +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(t_sel_e -
                                                                     t_sel_s)
                    .count();
        }

        // process collision (timed)
        std::tuple<bool, double3, double3> result;
        {
            auto t_proc_s = clock::now();
            result = process_collision(collision_type, is_electron, vcp, vn,
                                       mcp, mn);
            auto t_proc_e = clock::now();
            profiler.time_process_collision_ns +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(t_proc_e -
                                                                     t_proc_s)
                    .count();
        }

        auto t_end_total = clock::now();
        profiler.time_total_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_end_total -
                                                                 t_start_total)
                .count();
        return result;
    }

    // ----------------- PHYSICAL_ONLY or other physical scheme
    // ----------------- compute frequencies first (timed)
    double ion_freq = 0.0, cx_freq = 0.0;
    {
        auto t_cf0 = clock::now();
        std::tie(ion_freq, cx_freq) = compute_frequencies(vcp, vn, mcp, nn);
        auto t_cf1 = clock::now();
        profiler.time_compute_freq_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_cf1 - t_cf0)
                .count();
    }

    double total_freq = ion_freq + cx_freq;
    if (total_freq <= 0.0) {
        // record sample with P = 0 (frequencies were computed)
        profiler.add_freq_sample(ion_freq, cx_freq, 0.0);

        auto t_end_total = clock::now();
        profiler.time_total_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_end_total -
                                                                 t_start_total)
                .count();
        return {false, vcp, vn};
    }

    // compute collision probability (timed)
    double P_collision;
    {
        auto t_cp_s = clock::now();
        P_collision = compute_collision_probability(total_freq, dt);
        auto t_cp_e = clock::now();
        profiler.time_compute_prob_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_cp_e -
                                                                 t_cp_s)
                .count();
    }

    // record sample (we computed frequencies, use P_collision)
    profiler.add_freq_sample(ion_freq, cx_freq, P_collision);

    // check collision (timed)
    bool collided;
    {
        auto t_check_s = clock::now();
        collided = check_collision(P_collision);
        auto t_check_e = clock::now();
        profiler.time_check_collision_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_check_e -
                                                                 t_check_s)
                .count();
    }

    if (!collided) {
        auto t_end_total = clock::now();
        profiler.time_total_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_end_total -
                                                                 t_start_total)
                .count();
        return {false, vcp, vn};
    }

    // select collision type (timed)
    CollisionType collision_type;
    {
        auto t_sel_s = clock::now();
        collision_type =
            select_collision_type(is_electron, ion_freq, cx_freq, total_freq);
        auto t_sel_e = clock::now();
        profiler.time_select_type_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_sel_e -
                                                                 t_sel_s)
                .count();
    }

    // process collision (timed)
    std::tuple<bool, double3, double3> result;
    {
        auto t_proc_s = clock::now();
        result =
            process_collision(collision_type, is_electron, vcp, vn, mcp, mn);
        auto t_proc_e = clock::now();
        profiler.time_process_collision_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(t_proc_e -
                                                                 t_proc_s)
                .count();
    }

    auto t_end_total = clock::now();
    profiler.time_total_ns +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(t_end_total -
                                                             t_start_total)
            .count();

    return result;
}