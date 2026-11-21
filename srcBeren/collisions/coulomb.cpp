#include "collisions/coulomb.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

#include "core/vector3.h"
#include "runtime/env_options.h"
#include "runtime/logger.h"
#include "runtime/profiler.h"
#include "runtime/random.h"
#include "utils/sgs.h"

static inline double get_reduced_mass(double m1, double m2) {
    return m1 * m2 / (m1 + m2);
}

namespace {

// SplitMix64 finalizer: fast avalanche mixing of a 64-bit value.
inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// FNV-1a over a species name, with a separator so (a,b) and (ab,"") differ.
inline uint64_t mix_name(const std::string& name, uint64_t h) {
    for (unsigned char c : name) {
        h ^= c;
        h *= 0x00000100000001B3ULL;
    }
    h ^= 0x1F;
    h *= 0x00000100000001B3ULL;
    return h;
}

// Derive a stable per-channel seed from the operator's base seed and the
// ordered species names.  This is what makes e-e, i-i and e-i use different
// streams (previously all operators started from seed 42) without relying on
// construction order or thread count.
inline uint64_t derive_stream_seed(uint64_t baseSeed, const std::string& name1, const std::string& name2) {
    const uint64_t channelHash = mix_name(name2, mix_name(name1, 0xCBF29CE484222325ULL));
    return splitmix64(channelHash ^ splitmix64(baseSeed));
}

// Per-cell, per-apply seed.  Mixing the apply index keeps successive timesteps
// statistically independent while remaining deterministic for any thread count.
inline uint64_t cell_stream_seed(uint64_t opSeed, uint64_t applyIndex, int cellIndex) {
    uint64_t x = opSeed ^ splitmix64(applyIndex + 0x9E3779B97F4A7C15ULL);
    return splitmix64(x + 0xD1B54A32D192ED03ULL * static_cast<uint64_t>(cellIndex + 1));
}

// Read once; disabled statistics add only predictable branches to the pair loop.
const bool collisionStatsEnabled = envOptions::getenvParsed<bool>("COLLISION_STATS", false);

// Aggregate per-thread counters and print one compact line.  Thread times are
// recorded after each thread's `omp for nowait` chunk and before the implicit
// barrier at the end of the parallel region, so they exclude waiting (the flat
// timers, by contrast, include it).
void log_collision_stats(const std::string& label, const std::vector<CoulombCollisionStats>& threadStats,
                         const std::vector<double>& threadTimes, int nthreads, int numCells) {
    if (nthreads < 1)
        return;

    CoulombCollisionStats total;
    for (int t = 0; t < nthreads; ++t) {
        const CoulombCollisionStats& s = threadStats[t];
        total.pairs += s.pairs;
        total.variance_ge1 += s.variance_ge1;
        total.modu_small += s.modu_small;
        total.sum_var += s.sum_var;
        total.max_var = std::max(total.max_var, s.max_var);
        total.max_pairs_per_cell = std::max(total.max_pairs_per_cell, s.max_pairs_per_cell);
        total.cells_skipped += s.cells_skipped;
    }

    std::vector<double> times(threadTimes.begin(), threadTimes.begin() + nthreads);
    std::sort(times.begin(), times.end());
    const double tMin = times.front();
    const double tMax = times.back();
    const double tMedian = (nthreads % 2) ? times[nthreads / 2] : 0.5 * (times[nthreads / 2 - 1] + times[nthreads / 2]);

    std::ostringstream os;
    os << std::scientific << std::setprecision(3);
    os << "COLLISION_STATS " << label << " cells=" << numCells << " threads=" << nthreads << " pairs=" << total.pairs
       << " var>=1=" << total.variance_ge1 << " modu<1e-30=" << total.modu_small << " sum_var=" << total.sum_var
       << " max_var=" << total.max_var << " max_pairs_per_cell=" << total.max_pairs_per_cell
       << " cells_skipped=" << total.cells_skipped << " t_min=" << tMin << " t_median=" << tMedian << " t_max=" << tMax
       << " max/median=" << (tMedian > 0.0 ? tMax / tMedian : 0.0);
    logger::info(os.str());
}

// Momentum and kinetic energy of the particles of two species in one cell.
struct CellMoments {
    Vector3R momentum{0, 0, 0};
    double energy1 = 0;
    double energy2 = 0;
};

CellMoments cell_moments(const std::vector<Particle>& cell1, const std::vector<Particle>& cell2, double mass1,
                         double mass2) {
    CellMoments result;
    for (const Particle& p : cell1) {
        result.momentum += mass1 * p.velocity;
        result.energy1 += 0.5 * mass1 * p.velocity.dot(p.velocity);
    }
    for (const Particle& p : cell2) {
        result.momentum += mass2 * p.velocity;
        result.energy2 += 0.5 * mass2 * p.velocity.dot(p.velocity);
    }
    return result;
}

// Remove the energy error `excess` from one species by zero-angle inelastic scattering of random pairs:
// each pair keeps its momentum and changes its centre-of-mass energy by at most `fraction` of it.
// Returns the part of `excess` that could not be absorbed (fewer than two moving particles).
double absorb_energy(std::vector<Particle>& cell, double mass, double excess, double fraction, LehmerEngine& rndEng) {
    const int n = std::ssize(cell);
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    for (int pass = 0; pass < 1000 && excess != 0 && n >= 2; ++pass) {
        std::shuffle(order.begin(), order.end(), rndEng);
        bool changed = false;
        for (int k = 0; k + 1 < n && excess != 0; k += 2) {
            Particle& a = cell[order[k]];
            Particle& b = cell[order[k + 1]];
            const Vector3R relative = a.velocity - b.velocity;
            const double cm_energy = 0.25 * mass * relative.dot(relative);   // reduced mass m/2
            if (cm_energy == 0)
                continue;
            const double change = std::copysign(std::min(std::abs(excess), fraction * cm_energy), excess);
            const Vector3R shift = 0.5 * (std::sqrt(1.0 - change / cm_energy) - 1.0) * relative;
            a.velocity += shift;
            b.velocity -= shift;
            excess -= change;
            changed = true;
        }
        if (!changed)
            break;
    }
    return excess;
}

// Algorithm 1 of Angus et al. (2024), steps 7-14: shift velocities in proportion to the particle weight to
// restore the cell momentum, then share the energy error between the species by weight times energy and
// absorb it pairwise. Particle masses are species mass times weight.
void restore_moments(std::vector<Particle>& cell1, std::vector<Particle>& cell2, double mass1, double mass2, double w1,
                     double w2, const CellMoments& before, double fraction, LehmerEngine& rndEng) {
    const CellMoments scattered = cell_moments(cell1, cell2, mass1, mass2);
    const double weighted_mass = std::ssize(cell1) * mass1 * w1 + std::ssize(cell2) * mass2 * w2;
    const Vector3R shift = (scattered.momentum - before.momentum) / weighted_mass;
    for (Particle& p : cell1) p.velocity -= w1 * shift;
    for (Particle& p : cell2) p.velocity -= w2 * shift;
    const CellMoments corrected = cell_moments(cell1, cell2, mass1, mass2);
    const double excess = corrected.energy1 + corrected.energy2 - before.energy1 - before.energy2;
    const double share1 = w1 * corrected.energy1, share2 = w2 * corrected.energy2;
    if (excess == 0 || share1 + share2 == 0)
        return;
    const double excess1 = std::ssize(cell2) < 2   ? excess
                           : std::ssize(cell1) < 2 ? 0
                                                   : excess * share1 / (share1 + share2);
    double excess2 = excess - excess1;
    excess2 += absorb_energy(cell1, mass1, excess1, fraction, rndEng);
    absorb_energy(cell2, mass2, excess2, fraction, rndEng);   // a cell with one particle per species keeps the rest
}

}   // namespace

// ============================================================================
// CoulombCollisionOperator
// ============================================================================

CoulombCollisionOperator::CoulombCollisionOperator(const std::string& species1, const std::string& species2, double n0,
                                                   double coulomb_log, double energy_correction)
    : species1_name_(species1),
      species2_name_(species2),
      is_same_type_(species1 == species2),
      n0_(n0),
      coulomb_log_(coulomb_log),
      energy_correction_(energy_correction),
      stream_seed_(derive_stream_seed(envOptions::collisionSeed(), species1, species2)) {
    if (!std::isfinite(n0_) || n0_ <= 0 || !std::isfinite(coulomb_log_) || coulomb_log_ < 0)
        throw std::invalid_argument("Coulomb collision requires n0 > 0 and finite nonnegative lnL");
    if (!(energy_correction_ > 0 && energy_correction_ < 1))
        throw std::invalid_argument("Coulomb energy_correction_fraction must lie in (0, 1)");
}

std::string CoulombCollisionOperator::info() const {
    if (is_same_type_) {
        return "Coulomb(" + species1_name_ + "-" + species2_name_ + ", lnL=" + std::to_string(coulomb_log_) + ")";
    }
    return "Coulomb(" + species1_name_ + "-" + species2_name_ + ", lnL=" + std::to_string(coulomb_log_) + ")";
}

Vector3R CoulombCollisionOperator::scatter(const Vector3R& v1, const Vector3R& v2, double variance_numerator,
                                           double plasma_factor, double denominator_const, double variance_factor,
                                           std::uniform_real_distribution<double>& uniformDistr,
                                           std::normal_distribution<double>& normalDistr, LehmerEngine& rndEng,
                                           CoulombCollisionStats* stats) {
    const Vector3R u = v1 - v2;
    // In code units velocities are O(1), far from the range where hypot's
    // scaling matters; sqrt(dot) avoids hypot's divisions.
    const double modu = std::sqrt(u.dot(u));
    if (modu < 1.e-30) {
        if (stats) {
            ++stats->pairs;
            ++stats->modu_small;
        }
        return Vector3R(0, 0, 0);
    }

    // Preserve the arithmetic order while hoisting cell/operator constants:
    //   vf * (A * (numerator / (C * u^3)))
    // A, C, numerator and vf are hoisted to operator/cell/pair-class level.
    const double denominator = denominator_const * modu * modu * modu;
    const double variance = variance_factor * (plasma_factor * (variance_numerator / denominator));
    if (stats) {
        ++stats->pairs;
        stats->sum_var += variance;
        stats->max_var = std::max(stats->max_var, variance);
        if (!(variance < 1.0))
            ++stats->variance_ge1;
    }

    // Reuse one distribution object per thread.  param()/reset() reproduce the
    // state of a freshly constructed normal_distribution, so the RNG sequence
    // is identical to the old per-pair construction.
    normalDistr.param(std::normal_distribution<double>::param_type(0.0, std::sqrt(variance)));
    normalDistr.reset();
    const double sigma = (variance < 1.0) ? normalDistr(rndEng) : M_PI * uniformDistr(rndEng);
    const double sint = 2.0 * sigma / (1.0 + sigma * sigma);
    const double cost = 1.0 - 2.0 * sigma * sigma / (1.0 + sigma * sigma);

    const double phi = 2.0 * M_PI * uniformDistr(rndEng);
    const double cosp = cos(phi);
    const double sinp = sin(phi);
    const double up = sqrt(u.x() * u.x() + u.y() * u.y());
    const double one_minus_cost = 1.0 - cost;

    Vector3R du;
    // Two branches for the Takizuka-Abe rotation:
    //
    // Branch 1 (up ≈ 0): u is nearly parallel to z — rotate in XY plane.
    //   du_xy = modu * sin(θ) * (cos φ, sin φ),  du_z = -sign(u_z) * modu * (1 - cos θ)
    //   Result: |u'| = |u|, u' is rotated by θ isotropically around z.
    //
    // Branch 2 (general 3D): standard rotation of u by θ around an axis
    //   perpendicular to u at azimuth φ.  When up → 0, branch 2 approaches:
    //     dux ≈ modu * sin θ * cos(α + φ),  duy ≈ modu * sin θ * sin(α + φ)
    //   where α is the azimuth of the vanishing perpendicular component.
    //   Since φ ~ Uniform[0, 2π), the phase shift α is irrelevant —
    //   both branches give the same distribution.
    if (up < 1.e-16 * modu) {
        double sign_z = (u.z() >= 0.0) ? 1.0 : -1.0;
        du = Vector3R(modu * sint * cosp, modu * sint * sinp, -sign_z * modu * one_minus_cost);
    } else {
        // Hoist the two quotients (each is used by both components) instead of
        // doing four per-pair divisions.
        const double ux_over_up = u.x() / up;
        const double uy_over_up = u.y() / up;
        du = Vector3R(ux_over_up * u.z() * sint * cosp - uy_over_up * modu * sint * sinp - u.x() * one_minus_cost,
                      uy_over_up * u.z() * sint * cosp + ux_over_up * modu * sint * sinp - u.y() * one_minus_cost,
                      -up * sint * cosp - u.z() * one_minus_cost);
    }

    return du;
}

void CoulombCollisionOperator::collide_same_type(ParticlesArray& sp, double dt) {
    RECORD_TIMER;
    if (!std::isfinite(dt) || dt < 0)
        throw std::invalid_argument("Coulomb timestep must be finite and nonnegative");
    if (dt == 0)
        return;
    const double q = sp.charge();
    const double m1 = sp.mass();
    const int num_cells = sp.size();

    // Operator-level constants: previously recomputed (including pow and
    // several divisions) for every pair.
    const double wp3 = std::pow(SGS::get_plasma_freq(n0_), 3);
    constexpr double c3 = SGS::c * SGS::c * SGS::c;
    const double m = get_reduced_mass(m1, m1);
    const double plasma_factor = wp3 / (c3 * n0_);
    const double denominator_const = 8.0 * M_PI * m * m;
    const double charge_factor = coulomb_log_ * q * q * q * q;
    const double mu_over_m1 = m / m1;
    const double mu_over_m2 = m / m1;

    const bool statsOn = collisionStatsEnabled;
    const int maxThreads = statsOn ? omp_get_max_threads() : 0;
    std::vector<CoulombCollisionStats> threadStats(maxThreads);
    std::vector<double> threadTimes(maxThreads);
    int threadsUsed = 1;

#pragma omp parallel
    {
        timer::FlatTimer ompTimer("OMP section");

        const int tid = omp_get_thread_num();
        CoulombCollisionStats threadStat;
        CoulombCollisionStats* statsPtr = statsOn ? &threadStat : nullptr;
        const double tStart = statsOn ? omp_get_wtime() : 0.0;

        std::uniform_real_distribution<double> uniformDistr(0.0, 1.0);
        std::normal_distribution<double> normalDistr;
        std::vector<int> indices;

#pragma omp for schedule(static) nowait
        for (int pk = 0; pk < num_cells; pk++) {
            std::vector<Particle>& cell = sp(pk);
            const int N = cell.size();
            if (N < 2) {
                if (statsOn)
                    ++threadStat.cells_skipped;
                continue;
            }

            // Independent per-cell stream from a hash of (channel, apply, cell).
            LehmerEngine rndEng(cell_stream_seed(stream_seed_, apply_index_, pk));

            indices.resize(N);
            std::iota(indices.begin(), indices.end(), 0);

            const bool is_odd = (N % 2 == 1);
            const double n_density = N * sp.mpw();
            const double variance_numerator = charge_factor * n_density * dt;

            if (statsOn) {
                const int64_t cellPairs = is_odd ? int64_t{3} + (N - 3) / 2 : N / 2;
                threadStat.max_pairs_per_cell = std::max(threadStat.max_pairs_per_cell, cellPairs);
            }

            if (is_odd) {
                // Case 1b: first 3 particles form 3 pairs with variance_factor=0.5
                // Lazy Fisher-Yates for first 3 elements
                for (int k = 0; k < 3 && k < N; k++) {
                    std::uniform_int_distribution<int> dist(k, N - 1);
                    std::swap(indices[k], indices[dist(rndEng)]);
                }

                auto do_collide = [&](int i1, int i2, double variance_factor) {
                    const Vector3R du =
                        scatter(cell[i1].velocity, cell[i2].velocity, variance_numerator, plasma_factor,
                                denominator_const, variance_factor, uniformDistr, normalDistr, rndEng, statsPtr);
                    cell[i1].velocity += mu_over_m1 * du;
                    cell[i2].velocity -= mu_over_m2 * du;
                };

                do_collide(indices[0], indices[1], 0.5);
                do_collide(indices[1], indices[2], 0.5);
                do_collide(indices[0], indices[2], 0.5);

                // Remaining pairs: lazy Fisher-Yates from index 3 onward
                for (int k = 3; k + 1 < N; k += 2) {
                    std::uniform_int_distribution<int> dist1(k, N - 1);
                    std::swap(indices[k], indices[dist1(rndEng)]);
                    std::uniform_int_distribution<int> dist2(k + 1, N - 1);
                    std::swap(indices[k + 1], indices[dist2(rndEng)]);
                    do_collide(indices[k], indices[k + 1], 1.0);
                }
            } else {
                // Case 1a: all pairs with variance_factor=1.0
                // Lazy Fisher-Yates: shuffle only as needed
                for (int k = 0; k + 1 < N; k += 2) {
                    std::uniform_int_distribution<int> dist1(k, N - 1);
                    std::swap(indices[k], indices[dist1(rndEng)]);
                    std::uniform_int_distribution<int> dist2(k + 1, N - 1);
                    std::swap(indices[k + 1], indices[dist2(rndEng)]);

                    const Vector3R du =
                        scatter(cell[indices[k]].velocity, cell[indices[k + 1]].velocity, variance_numerator,
                                plasma_factor, denominator_const, 1.0, uniformDistr, normalDistr, rndEng, statsPtr);
                    cell[indices[k]].velocity += mu_over_m1 * du;
                    cell[indices[k + 1]].velocity -= mu_over_m2 * du;
                }
            }
        }

        // Record this thread's chunk time before the end-of-parallel barrier:
        // with nowait this excludes the time spent waiting for the slowest thread.
        if (statsOn) {
            threadTimes[tid] = omp_get_wtime() - tStart;
            threadStats[tid] = threadStat;
#pragma omp master
            threadsUsed = omp_get_num_threads();
        }
    }

    if (statsOn)
        log_collision_stats(info(), threadStats, threadTimes, threadsUsed, num_cells);

    // No master advance: the per-cell seed already depends on apply_index_.
    ++apply_index_;
}

void CoulombCollisionOperator::collide_diff_type(ParticlesArray& sp1, ParticlesArray& sp2, double dt) {
    RECORD_TIMER;
    if (!std::isfinite(dt) || dt < 0)
        throw std::invalid_argument("Coulomb timestep must be finite and nonnegative");
    if (sp1.cells().size() != sp2.cells().size())
        throw std::invalid_argument("Coulomb species must use the same grid");
    // Unequal weights (Angus et al. 2024): the density uses the larger weight, the heavier-weight particle
    // scatters with probability w_min / w_max, and each cell then restores momentum and energy exactly.
    const double w1 = sp1.mpw();
    const double w2 = sp2.mpw();
    const bool weighted = w1 != w2;
    const double w_max = std::max(w1, w2);
    const double update1 = w1 > w2 ? w2 / w1 : 1.0;
    const double update2 = w2 > w1 ? w1 / w2 : 1.0;
    if (dt == 0)
        return;
    const double q1 = sp1.charge();
    const double m1 = sp1.mass();
    const double q2 = sp2.charge();
    const double m2 = sp2.mass();
    const int num_cells = sp1.size();

    // Operator-level constants: previously recomputed (including pow and
    // several divisions) for every pair.
    const double wp3 = std::pow(SGS::get_plasma_freq(n0_), 3);
    constexpr double c3 = SGS::c * SGS::c * SGS::c;
    const double m = get_reduced_mass(m1, m2);
    const double plasma_factor = wp3 / (c3 * n0_);
    const double denominator_const = 8.0 * M_PI * m * m;
    const double charge_factor = coulomb_log_ * q1 * q1 * q2 * q2;
    const double mu_over_m1 = m / m1;
    const double mu_over_m2 = m / m2;

    const bool statsOn = collisionStatsEnabled;
    const int maxThreads = statsOn ? omp_get_max_threads() : 0;
    std::vector<CoulombCollisionStats> threadStats(maxThreads);
    std::vector<double> threadTimes(maxThreads);
    int threadsUsed = 1;

#pragma omp parallel
    {
        timer::FlatTimer ompTimer("OMP section");

        const int tid = omp_get_thread_num();
        CoulombCollisionStats threadStat;
        CoulombCollisionStats* statsPtr = statsOn ? &threadStat : nullptr;
        const double tStart = statsOn ? omp_get_wtime() : 0.0;

        std::uniform_real_distribution<double> uniformDistr(0.0, 1.0);
        std::normal_distribution<double> normalDistr;
        std::vector<int> indices_large;
        std::vector<int> indices_small;

#pragma omp for schedule(static) nowait
        for (int pk = 0; pk < num_cells; pk++) {
            std::vector<Particle>& cell1 = sp1(pk);
            std::vector<Particle>& cell2 = sp2(pk);
            const int N1 = cell1.size();
            const int N2 = cell2.size();
            if (N1 == 0 || N2 == 0) {
                if (statsOn)
                    ++threadStat.cells_skipped;
                continue;
            }

            // Independent per-cell stream from a hash of (channel, apply, cell).
            LehmerEngine rndEng(cell_stream_seed(stream_seed_, apply_index_, pk));

            const bool sp1_is_larger = (N1 >= N2);
            const int N_large = sp1_is_larger ? N1 : N2;
            const int N_small = sp1_is_larger ? N2 : N1;

            if (statsOn)
                threadStat.max_pairs_per_cell = std::max(threadStat.max_pairs_per_cell, int64_t{N_large});

            indices_large.resize(N_large);
            std::iota(indices_large.begin(), indices_large.end(), 0);

            indices_small.resize(N_small);
            std::iota(indices_small.begin(), indices_small.end(), 0);

            // Lazy Fisher-Yates for the large array
            // For the small array, we cycle through deterministically (Takizuka Case 2)
            const int r = N_large % N_small;
            const int quot = N_large / N_small;
            const int first_group_size = r * (quot + 1);

            // Shuffle the large array lazily as we iterate
            // Shuffle the small array once (full, since all get reused)
            for (int k = 0; k < N_small; k++) {
                std::uniform_int_distribution<int> dist(k, N_small - 1);
                std::swap(indices_small[k], indices_small[dist(rndEng)]);
            }

            // Equals min(N1 w1, N2 w2) bit for bit when the weights are equal.
            const double variance_numerator = charge_factor * (w_max * std::min(N1, N2)) * dt;
            const CellMoments before = weighted ? cell_moments(cell1, cell2, m1 * w1, m2 * w2) : CellMoments{};

            for (int first_ind = 0; first_ind < N_large; first_ind++) {
                // Lazy Fisher-Yates for large array
                std::uniform_int_distribution<int> dist(first_ind, N_large - 1);
                std::swap(indices_large[first_ind], indices_large[dist(rndEng)]);

                // Determine partner from small array (Takizuka Case 2 mapping)
                int second_ind;
                if (first_ind < first_group_size) {
                    second_ind = first_ind / (quot + 1);
                } else {
                    second_ind = r + (first_ind - first_group_size) / quot;
                }

                int idx1, idx2;
                if (sp1_is_larger) {
                    idx1 = indices_large[first_ind];
                    idx2 = indices_small[second_ind];
                } else {
                    idx1 = indices_small[second_ind];
                    idx2 = indices_large[first_ind];
                }

                const Vector3R du =
                    scatter(cell1[idx1].velocity, cell2[idx2].velocity, variance_numerator, plasma_factor,
                            denominator_const, 1.0, uniformDistr, normalDistr, rndEng, statsPtr);
                if (!weighted) {
                    cell1[idx1].velocity += mu_over_m1 * du;
                    cell2[idx2].velocity -= mu_over_m2 * du;
                    continue;
                }
                const double draw = uniformDistr(rndEng);
                if (update1 == 1.0 || draw < update1)
                    cell1[idx1].velocity += mu_over_m1 * du;
                if (update2 == 1.0 || draw < update2)
                    cell2[idx2].velocity -= mu_over_m2 * du;
            }
            if (weighted)
                restore_moments(cell1, cell2, m1 * w1, m2 * w2, w1, w2, before, energy_correction_, rndEng);
        }

        if (statsOn) {
            threadTimes[tid] = omp_get_wtime() - tStart;
            threadStats[tid] = threadStat;
#pragma omp master
            threadsUsed = omp_get_num_threads();
        }
    }

    if (statsOn)
        log_collision_stats(info(), threadStats, threadTimes, threadsUsed, num_cells);

    // No master advance: the per-cell seed already depends on apply_index_.
    ++apply_index_;
}

void CoulombCollisionOperator::apply(Species& species, double dt) {
    RECORD_TIMER;
    if (is_same_type_) {
        ParticlesArray* sp = find_species(species, species1_name_);
        if (!sp) {
            std::cerr << "CoulombCollision: species '" << species1_name_ << "' not found\n";
            return;
        }
        collide_same_type(*sp, dt);
    } else {
        ParticlesArray* sp1 = find_species(species, species1_name_);
        ParticlesArray* sp2 = find_species(species, species2_name_);
        if (!sp1 || !sp2) {
            std::cerr << "CoulombCollision: species '" << species1_name_ << "' or '" << species2_name_
                      << "' not found\n";
            return;
        }
        collide_diff_type(*sp1, *sp2, dt);
    }
}
