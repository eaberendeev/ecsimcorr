#pragma once

#include <algorithm>
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

#include "cross_section.h"
#include "sgs.h"
#include "vector3.h"

enum class CollisionScheme { NULL_COLLISION, PHYSICAL_ONLY };

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
    double sum_ionization_freq = 0.0;
    double sum_cx_freq = 0.0;
    double sum_P_collision = 0.0;
    double max_ionization_freq = 0.0;
    double max_cx_freq = 0.0;
    double max_P_collision = 0.0;

    double sum_ionization_sigma = 0.0;
    double sum_cx_sigma = 0.0;
    double max_ionization_sigma = 0.0;
    double max_cx_sigma = 0.0;

    // counters for non-zero frequencies
    std::uint64_t count_nonzero_ionization = 0;
    std::uint64_t count_nonzero_ionization_sigma = 0;
    std::uint64_t count_nonzero_cx = 0;
    std::uint64_t count_nonzero_P = 0;

    void add_freq_sample(double ionization_freq, double cx_freq, double P_collision) {
        freq_samples += 1;
        sum_ionization_freq += ionization_freq;
        sum_cx_freq += cx_freq;
        sum_P_collision += P_collision;
        max_ionization_freq = std::max(max_ionization_freq, ionization_freq);
        max_cx_freq = std::max(max_cx_freq, cx_freq);
        max_P_collision = std::max(max_P_collision, P_collision);
        count_nonzero_ionization += (ionization_freq > 0.0);
        count_nonzero_cx += (cx_freq > 0.0);
        count_nonzero_P += (P_collision > 0.0);
    }
    void add_sigma_sample(double ionization_sigma, double cx_sigma) {
        sum_ionization_sigma += ionization_sigma;
        sum_cx_sigma += cx_sigma;
        max_ionization_sigma = std::max(max_ionization_sigma, ionization_sigma);
        max_cx_sigma = std::max(max_cx_sigma, cx_sigma);
        count_nonzero_ionization_sigma += (ionization_sigma > 0.0);
    }
    void reset() {
        calls = freq_samples = 0;
        time_total_ns = time_compute_freq_ns = time_compute_prob_ns = time_check_collision_ns = 0;
        time_select_type_ns = time_process_collision_ns = 0;
        sum_ionization_freq = sum_cx_freq = sum_P_collision = 0.0;
        max_ionization_freq = max_cx_freq = max_P_collision = 0.0;
        sum_ionization_sigma = sum_cx_sigma = 0.0;
        max_ionization_sigma = max_cx_sigma = 0.0;
        count_nonzero_ionization = count_nonzero_ionization_sigma = count_nonzero_cx = count_nonzero_P = 0;
    }

    void print_report(std::ostream& os = std::cout) const {
        os << std::scientific << std::setprecision(12);
        os << "Collision profiler report:\n";
        os << "  total calls: " << calls << "\n";
        os << "  freq_samples (compute_frequencies called): " << freq_samples << "\n";
        auto ns_to_ms = [](std::uint64_t ns) { return double(ns) / 1e6; };
        os << "  total time (ms)           : " << ns_to_ms(time_total_ns) << "\n";
        os << "  compute_frequencies (ms)  : " << ns_to_ms(time_compute_freq_ns) << "\n";
        os << "  compute_collision_prob (ms): " << ns_to_ms(time_compute_prob_ns) << "\n";
        os << "  check_collision (ms)      : " << ns_to_ms(time_check_collision_ns) << "\n";
        os << "  select_type (ms)          : " << ns_to_ms(time_select_type_ns) << "\n";
        os << "  process_collision (ms)    : " << ns_to_ms(time_process_collision_ns) << "\n";
        if (freq_samples > 0) {
            os << "  avg ionization_freq : " << (sum_ionization_freq / double(freq_samples))
               << ", max: " << max_ionization_freq << ", nonzero count: " << count_nonzero_ionization << "\n";
            os << "  avg cx_freq  : " << (sum_cx_freq / double(freq_samples)) << ", max: " << max_cx_freq
               << ", nonzero count: " << count_nonzero_cx << "\n";
            os << "  avg ionization_sigma : " << (sum_ionization_sigma / double(freq_samples))
               << ", max: " << max_ionization_sigma << ", nonzero count: " << count_nonzero_ionization_sigma << "\n";
            os << "  avg cx_sigma  : " << (sum_cx_sigma / double(freq_samples)) << ", max: " << max_cx_sigma
               << ", nonzero count: " << count_nonzero_cx << "\n";
            os << "  avg P_coll   : " << (sum_P_collision / double(freq_samples)) << ", max: " << max_P_collision
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
    explicit ColliderWithNeutrals(const double _n0, CollisionScheme scheme = CollisionScheme::PHYSICAL_ONLY,
                                  CollisionProcessOptions process_opts = CollisionProcessOptions()) {
        n0 = _n0;
        wp = SGS::get_plasma_freq(_n0);
        l0 = SGS::c / wp;
        dpd0 = n0 * l0;
        scheme_mode = scheme;
        process_options = process_opts;
        uint64_t seed = static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
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

    void set_scheme(CollisionScheme new_scheme) {
        scheme_mode = new_scheme;
    }
    CollisionScheme get_scheme() const {
        return scheme_mode;
    }
    void set_process_options(const CollisionProcessOptions& new_options) {
        process_options = new_options;
    }
    const CollisionProcessOptions& get_process_options() const {
        return process_options;
    }
    std::tuple<bool, Vector3R, Vector3R> collision_with_neutral(Vector3R& vcp, Vector3R& vn, double mcp, double mn,
                                                                double nn, double dt, double freq_max);

    double total_collision_frequency(const Vector3R& vcp, const Vector3R& vn, double mcp, double nn);

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
        uint64_t r = rng();                                       // быстрый вызов генератора
        uint64_t top53 = r >> 11;                                 // получаем 53 бита
        constexpr double inv_2pow53 = 1.0 / 9007199254740992.0;   // 1/2^53
        return double(top53) * inv_2pow53;
    }
    CollisionProcessOptions process_options;

    std::pair<double, double> compute_frequencies(const Vector3R& vcp, const Vector3R& vn, double mcp, double nn);

    bool check_collision(double P_collision);

    CollisionType select_collision_type(bool is_electron, double ionization_freq, double cx_freq, double freq_bound);

    CollisionScheme scheme_mode;
};
