/**
 * @file collisions_with_neutrals.h
 * @author Morozov O. P.
 * @brief Declares the neutral collision handler and related configuration
 * types.
 */
#pragma once

#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>
#include <utility>

#include "Vec.h"
#include "cross_section.h"
#include "sgs.h"

enum class CollisionScheme {
    NULL_COLLISION,
    PHYSICAL_ONLY
};

struct CollisionProcessOptions {
    bool electron_ionization = true;
    bool proton_ionization = true;
    bool proton_charge_exchange = true;
};

struct CollisionProfiler {
    std::uint64_t calls = 0;
    std::uint64_t freq_samples = 0;

    // times (ns)
    std::uint64_t time_total_ns = 0;
    std::uint64_t time_compute_freq_ns = 0;
    std::uint64_t time_compute_prob_ns = 0;
    std::uint64_t time_check_collision_ns = 0;
    std::uint64_t time_select_type_ns = 0;
    std::uint64_t time_process_collision_ns = 0;

    // sum and max values (by freq_samples)
    double sum_ion_freq = 0.0;
    double sum_cx_freq = 0.0;
    double sum_P_collision = 0.0;
    double max_ion_freq = 0.0;
    double max_cx_freq = 0.0;
    double max_P_collision = 0.0;

    // counters for non-zero frequencies
    std::uint64_t count_nonzero_ion = 0;
    std::uint64_t count_nonzero_cx = 0;
    std::uint64_t count_nonzero_P = 0;

    void add_freq_sample(double ion, double cx, double P_collision) {
        freq_samples += 1;
        sum_ion_freq += ion;
        sum_cx_freq += cx;
        sum_P_collision += P_collision;
        if (ion > max_ion_freq)
            max_ion_freq = ion;
        if (cx > max_cx_freq)
            max_cx_freq = cx;
        if (P_collision > max_P_collision)
            max_P_collision = P_collision;
        if (ion > 0.0)
            ++count_nonzero_ion;
        if (cx > 0.0)
            ++count_nonzero_cx;
        if (P_collision > 0.0)
            ++count_nonzero_P;
    }

    void reset() {
        calls = freq_samples = 0;
        time_total_ns = time_compute_freq_ns = time_compute_prob_ns =
            time_check_collision_ns = 0;
        time_select_type_ns = time_process_collision_ns = 0;
        sum_ion_freq = sum_cx_freq = sum_P_collision = 0.0;
        max_ion_freq = max_cx_freq = max_P_collision = 0.0;
        count_nonzero_ion = count_nonzero_cx = count_nonzero_P = 0;
    }

    void print_report(std::ostream& os = std::cout) const {
        os << std::scientific << std::setprecision(12);
        os << "Collision profiler report:\n";
        os << "  total calls: " << calls << "\n";
        os << "  freq_samples (compute_frequencies called): " << freq_samples
           << "\n";
        auto ns_to_ms = [](std::uint64_t ns) { return double(ns) / 1e6; };
        os << "  total time (ms)           : " << ns_to_ms(time_total_ns)
           << "\n";
        os << "  compute_frequencies (ms)  : " << ns_to_ms(time_compute_freq_ns)
           << "\n";
        os << "  compute_collision_prob (ms): "
           << ns_to_ms(time_compute_prob_ns) << "\n";
        os << "  check_collision (ms)      : "
           << ns_to_ms(time_check_collision_ns) << "\n";
        os << "  select_type (ms)          : " << ns_to_ms(time_select_type_ns)
           << "\n";
        os << "  process_collision (ms)    : "
           << ns_to_ms(time_process_collision_ns) << "\n";
        if (freq_samples > 0) {
            os << "  avg ion_freq : " << (sum_ion_freq / double(freq_samples))
               << ", max: " << max_ion_freq
               << ", nonzero count: " << count_nonzero_ion << "\n";
            os << "  avg cx_freq  : " << (sum_cx_freq / double(freq_samples))
               << ", max: " << max_cx_freq
               << ", nonzero count: " << count_nonzero_cx << "\n";
            os << "  avg P_coll   : "
               << (sum_P_collision / double(freq_samples))
               << ", max: " << max_P_collision
               << ", nonzero count: " << count_nonzero_P << "\n";
            double frac = double(freq_samples) / double(calls);
            os << "  fraction freq_samples/calls: " << frac << "\n";
        } else {
            os << "  (no frequency samples collected)\n";
        }
    }
};

class ColliderWithNeutrals {
   public:
    explicit ColliderWithNeutrals(
        const double _n0,
        CollisionScheme scheme = CollisionScheme::PHYSICAL_ONLY,
        CollisionProcessOptions process_opts = CollisionProcessOptions()) {
        n0 = _n0;
        wp = SGS::get_plasma_freq(_n0);
        l0 = SGS::c / wp;
        dpd0 = n0 * l0;
        scheme_mode = scheme;
        process_options = process_opts;
        uint64_t seed =
            static_cast<uint64_t>(std::chrono::high_resolution_clock::now()
                                      .time_since_epoch()
                                      .count());
        rng.seed(seed);
#pragma omp master
        {
            std::cout << "ColliderWithNeutrals create with parametrs: \n";
            std::cout << "n0: " << _n0 << "\n";
            std::cout << "wp: " << wp << "\n";
            std::cout << "l0: " << l0 << "\n";
            std::cout << "dpd0: " << dpd0 << "\n";
        }
    }

    void set_scheme(CollisionScheme new_scheme) { scheme_mode = new_scheme; }
    CollisionScheme get_scheme() const { return scheme_mode; }
    void set_process_options(const CollisionProcessOptions& new_options) {
        process_options = new_options;
    }
    const CollisionProcessOptions& get_process_options() const {
        return process_options;
    }
    std::tuple<bool, double3, double3> collision_with_neutral(
        double3& vcp, double3& vn, double mcp, double mn, double ncp, double nn,
        double dt, double freq_max);

    double total_collision_frequency(const double3& vcp, const double3& vn,
                                     double mcp, double nn) const;

    double Sigma_e(double E) const;
    double Sigma_p(double E) const;
    double Sigma_cx(double E) const;
    CollisionProfiler profiler;
    std::mt19937_64 rng;

   private:
    double n0;
    double wp;
    double l0;
    double dpd0;
    inline double rng_uniform01() {
        // Берём 53 старших бит из 64-битного случайного числа -> [0, 1)
        // Это стандартная техника для равномерного double с 53-битной
        // мантиссой.
        uint64_t r = rng();         // быстрый вызов генератора
        uint64_t top53 = r >> 11;   // получаем 53 бита
        constexpr double inv_2pow53 = 1.0 / 9007199254740992.0;   // 1/2^53
        return double(top53) * inv_2pow53;
    }
    CollisionProcessOptions process_options;

    std::pair<double, double> compute_frequencies(const double3& vcp,
                                                  const double3& vn, double mcp,
                                                  double nn) const;

    bool check_collision(double P_collision);

    CollisionType select_collision_type(bool is_electron, double ion_freq,
                                        double cx_freq, double freq_bound);

    CollisionScheme scheme_mode;
};
