#include "collision_manager.h"

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

#include "collisions_with_neutrals.h"
#include "core/Vector3.h"
#include "env_options.h"
#include "logger.h"
#include "random.h"
#include "timer.h"

// Jump-ahead stride separating per-cell random streams (see collision.cpp).
static constexpr int64_t collisionStreamStride = 1 << 16;

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

}   // namespace

// ============================================================================
// CoulombCollisionOperator
// ============================================================================

CoulombCollisionOperator::CoulombCollisionOperator(const std::string& species1, const std::string& species2, double n0,
                                                   double coulomb_log)
    : species1_name_(species1),
      species2_name_(species2),
      is_same_type_(species1 == species2),
      n0_(n0),
      coulomb_log_(coulomb_log),
      stream_seed_(derive_stream_seed(envOptions::collisionSeed(), species1, species2)) {
    if (!std::isfinite(n0_) || n0_ <= 0 || !std::isfinite(coulomb_log_) || coulomb_log_ < 0)
        throw std::invalid_argument("Coulomb collision requires n0 > 0 and finite nonnegative lnL");
}

std::string CoulombCollisionOperator::info() const {
    if (is_same_type_) {
        return "Coulomb(" + species1_name_ + "-" + species2_name_ + ", lnL=" + std::to_string(coulomb_log_) + ")";
    }
    return "Coulomb(" + species1_name_ + "-" + species2_name_ + ", lnL=" + std::to_string(coulomb_log_) + ")";
}

void CoulombCollisionOperator::bin_collide(Vector3R& v1, Vector3R& v2, double variance_numerator, double plasma_factor,
                                           double denominator_const, double variance_factor, double mu_over_m1,
                                           double mu_over_m2, std::uniform_real_distribution<double>& uniformDistr,
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
        return;
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

    v1 += mu_over_m1 * du;
    v2 -= mu_over_m2 * du;
}

void CoulombCollisionOperator::collide_same_type(ParticlesArray& sp, double dt) {
    RECORD_TIMER;
    if (!std::isfinite(dt) || dt < 0)
        throw std::invalid_argument("Coulomb timestep must be finite and nonnegative");
    if (dt == 0)
        return;
    const double q = sp.charge;
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
        timer::flatTimer ompTimer("OMP section");

        const int tid = omp_get_thread_num();
        CoulombCollisionStats threadStat;
        CoulombCollisionStats* statsPtr = statsOn ? &threadStat : nullptr;
        const double tStart = statsOn ? omp_get_wtime() : 0.0;

        std::uniform_real_distribution<double> uniformDistr(0.0, 1.0);
        std::normal_distribution<double> normalDistr;
        std::vector<int> indices;

#pragma omp for schedule(static) nowait
        for (int pk = 0; pk < num_cells; pk++) {
            auto& cell = sp.particlesData(pk);
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
                    bin_collide(cell[i1].velocity, cell[i2].velocity, variance_numerator, plasma_factor,
                                denominator_const, variance_factor, mu_over_m1, mu_over_m2, uniformDistr, normalDistr,
                                rndEng, statsPtr);
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

                    bin_collide(cell[indices[k]].velocity, cell[indices[k + 1]].velocity, variance_numerator,
                                plasma_factor, denominator_const, 1.0, mu_over_m1, mu_over_m2, uniformDistr,
                                normalDistr, rndEng, statsPtr);
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
    if (sp1.particlesData.size() != sp2.particlesData.size())
        throw std::invalid_argument("Coulomb species must use the same grid");
    // This binary rotation conserves physical momentum only for equal macro weights.
    const double weight_scale = std::max(std::abs(sp1.mpw()), std::abs(sp2.mpw()));
    if (std::abs(sp1.mpw() - sp2.mpw()) > 8 * std::numeric_limits<double>::epsilon() * weight_scale)
        throw std::invalid_argument("Coulomb interspecies collisions require equal macroparticle weights");
    if (dt == 0)
        return;
    const double q1 = sp1.charge;
    const double m1 = sp1.mass();
    const double q2 = sp2.charge;
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
        timer::flatTimer ompTimer("OMP section");

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
            auto& cell1 = sp1.particlesData(pk);
            auto& cell2 = sp2.particlesData(pk);
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

            const double n1_density = N1 * sp1.mpw();
            const double n2_density = N2 * sp2.mpw();
            const double variance_numerator = charge_factor * std::min(n1_density, n2_density) * dt;

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

                bin_collide(cell1[idx1].velocity, cell2[idx2].velocity, variance_numerator, plasma_factor,
                            denominator_const, 1.0, mu_over_m1, mu_over_m2, uniformDistr, normalDistr, rndEng,
                            statsPtr);
            }
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

void CoulombCollisionOperator::apply(Species& species, [[maybe_unused]] const Domain& domain, double dt) {
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

// ============================================================================
// NeutralCollisionOperator
// ============================================================================

NeutralCollisionOperator::NeutralCollisionOperator(const NeutralCollisionConfig& config, double n0)
    : config_(config), n0_(n0), baseRndEng_(77) {
}

std::string NeutralCollisionOperator::info() const {
    return "Neutral(" + config_.charged_species + "-" + config_.neutral_species + ")";
}

void NeutralCollisionOperator::apply(Species& species, const Domain& domain, double dt) {
    RECORD_TIMER;

    ParticlesArray* charged = find_species(species, config_.charged_species);
    ParticlesArray* neutrals = find_species(species, config_.neutral_species);
    ParticlesArray* electrons = find_species(species, config_.electron_product);
    ParticlesArray* ions = find_species(species, config_.ion_product);

    if (!charged || !neutrals || !electrons || !ions) {
        std::cerr << "NeutralCollision: missing species\n";
        return;
    }

    CollisionScheme scheme = CollisionScheme::PHYSICAL_ONLY;
    if (config_.scheme == "null_collision") {
        scheme = CollisionScheme::NULL_COLLISION;
    }

    CollisionProcessOptions opts;
    opts.electron_ionization = config_.electron_ionization;
    opts.proton_ionization = config_.proton_ionization;
    opts.proton_charge_exchange = config_.proton_charge_exchange;

    const double m1 = charged->mass();
    const double m2 = neutrals->mass();

#pragma omp parallel
    {
        timer::flatTimer ompTimer("OMP section");
        ColliderWithNeutrals collider(n0_, scheme, opts);

#pragma omp for schedule(static)
        for (int pk = 0; pk < charged->size(); pk++) {
            int pInCell = charged->particlesData(pk).size();
            int nInCell = neutrals->particlesData(pk).size();
            if (pInCell == 0 || nInCell == 0)
                continue;

            LehmerEngine rndEng = baseRndEng_;
            rndEng.discard(static_cast<int64_t>(pk) * collisionStreamStride);

            double n_neutrals = nInCell / (double) neutrals->NumPartPerCell;
            auto& neutrals_data = neutrals->particlesData(pk);
            int current_neutral_count = nInCell;

            for (int i = 0; i < pInCell && current_neutral_count > 0; i++) {
                std::uniform_int_distribution<> dis(0, current_neutral_count - 1);
                int randomIndex = dis(rndEng);

                Vector3R v1 = charged->particlesData(pk)[i].velocity;
                Vector3R v2 = neutrals_data[randomIndex].velocity;

                auto [is_collided, ve, vi] = collider.collision_with_neutral(v1, v2, m1, m2, n_neutrals, dt, 1.0 / dt);

                if (is_collided) {
                    Vector3R coord = neutrals_data[randomIndex].coord;
                    Particle pe(coord, ve);
                    Particle pi(coord, vi);

                    electrons->add_particle(pe);
                    ions->add_particle(pi);

                    charged->particlesData(pk)[i].velocity = v1;
                    std::swap(neutrals_data[randomIndex], neutrals_data[current_neutral_count - 1]);
                    current_neutral_count--;
                } else {
                    charged->particlesData(pk)[i].velocity = v1;
                    neutrals_data[randomIndex].velocity = v2;
                }
            }

            if (current_neutral_count < static_cast<int>(neutrals_data.size())) {
                neutrals_data.resize(current_neutral_count);
            }
        }
    }

    baseRndEng_.discard(charged->size() * collisionStreamStride);

    neutrals->move(dt);
    neutrals->update_cells(domain);
}

// ============================================================================
// CollisionManager
// ============================================================================

void CollisionManager::init_from_json(const nlohmann::json& config, double n0) {
    operators_.clear();

    if (!config.contains("Collisions")) {
        return;
    }

    const auto& collisions = config["Collisions"];
    if (!collisions.is_array()) {
        std::cerr << "CollisionManager: 'Collisions' must be an array\n";
        return;
    }

    for (const auto& entry : collisions) {
        if (!entry.contains("type")) {
            std::cerr << "CollisionManager: entry missing 'type'\n";
            continue;
        }

        const std::string type = entry["type"].get<std::string>();

        if (type == "coulomb") {
            std::string sp1 = entry.value("species1", "");
            std::string sp2 = entry.value("species2", "");
            double lnL = entry.value("coulomb_log", 15.0);
            if (sp1.empty() || sp2.empty()) {
                std::cerr << "CollisionManager: coulomb entry needs species1, species2\n";
                continue;
            }
            operators_.push_back(std::make_unique<CoulombCollisionOperator>(sp1, sp2, n0, lnL));
        } else if (type == "neutral_ionization") {
            NeutralCollisionConfig nc;
            nc.charged_species = entry.value("charged_species", "");
            nc.neutral_species = entry.value("neutral_species", "Neutrals");
            nc.electron_product = entry.value("electron_product", "Electrons");
            nc.ion_product = entry.value("ion_product", "Ions");
            nc.scheme = entry.value("scheme", "physical_only");
            nc.electron_ionization = entry.value("electron_ionization", true);
            nc.proton_ionization = entry.value("proton_ionization", true);
            nc.proton_charge_exchange = entry.value("proton_charge_exchange", true);
            if (nc.charged_species.empty()) {
                std::cerr << "CollisionManager: neutral_ionization needs charged_species\n";
                continue;
            }
            operators_.push_back(std::make_unique<NeutralCollisionOperator>(nc, n0));
        } else {
            std::cerr << "CollisionManager: unknown collision type '" << type << "'\n";
        }
    }

    if (operators_.empty())
        return;
    std::cout << "Collisions:\n";
    for (const auto& op : operators_) {
        std::cout << "  " << op->info() << "\n";
    }
}

void CollisionManager::apply(Species& species, const Domain& domain, double dt) {
    for (auto& op : operators_) {
        op->apply(species, domain, dt);
    }
}
