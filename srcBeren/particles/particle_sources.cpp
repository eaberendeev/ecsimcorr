#include "particle_sources.h"

#include <omp.h>

#include <atomic>
#include <iterator>

#include "particles/particle_array.h"
#include "runtime/profiler.h"
#include "utils/util.h"

InjectedParticles add_particles_from_distribution(ParticlesArray& species, IDistribution& dist, LehmerEngine& rng_space,
                                                  LehmerEngine& rng_momentum, const Domain& domain, double dt,
                                                  bool check_boundaries) {
    const int count = dist.get_count_to_inject();
    RECORD_TIMER_PARAMS(count);

    // spawn more particles than this in parallel; below that a serial loop is used
    constexpr int parallelSpawnThreshold = 15'000;
    // upper bound of RNG draws per particle (box-muller + rejection); used to keep per-particle streams disjoint
    constexpr int64_t largeNumber = 32;

    const bool useThreads = (count > parallelSpawnThreshold);

    // small fixed table of cell locks (hash of the cell index): avoids allocating
    // O(domain.total_size()) atomics on every call; collisions only serialize two cells
    constexpr int64_t locksCount = 1 << 16;
    std::vector<std::atomic<int64_t>> locks(useThreads ? locksCount : 0);
    if (useThreads) {
#pragma omp parallel for
        for (int i = 0; i < std::ssize(locks); ++i) {
            locks[i].store(0);
        }
    }

    // run-to run reproducible energy result
    const int nthr = useThreads ? omp_get_max_threads() : 1;
    std::vector<double> energyArray(nthr, 0.0);
    int accepted_count = 0;
#pragma omp parallel if (useThreads) num_threads(nthr) reduction(+ : accepted_count)
    {
        timer::flatTimer timer("OMP section");

        double threadLocalEnergy = 0.0;

#pragma omp for
        for (int i = 0; i < count; ++i) {
            LehmerEngine engSpace = rng_space;
            LehmerEngine engMomentum = rng_momentum;

            engSpace.discard(i * largeNumber);
            engMomentum.discard(i * largeNumber);

            Vector3R position = dist.sample_position(engSpace);
            const Vector3R velocity = dist.sample_velocity(position, engMomentum);

            if (dist.is_bound_injection()) {
                position += velocity * dt;
            }

            Particle particle(position, velocity);

            if (!check_boundaries || domain.contains(position)) {
                const Vector3I cell_id = species.get_domain().storage_cell_index(particle.coord);
                assert(species.get_domain().grid.contains_cell_index(cell_id));
                ++accepted_count;
                threadLocalEnergy += dist.get_energy(velocity);

                if (useThreads) {
                    std::atomic<int64_t>& lock =
                        locks[species.get_domain().sind(cell_id.x(), cell_id.y(), cell_id.z()) & (locksCount - 1)];
                    int64_t expected = 0;
                    while (!lock.compare_exchange_strong(expected, 1)) {
                        expected = 0;
                        asm volatile("pause");
                    }
                    species(cell_id.x(), cell_id.y(), cell_id.z()).push_back(particle);
                    lock.store(0);
                } else {
                    species.add_particle(particle);
                }
            }
        }

        energyArray[omp_get_thread_num()] = threadLocalEnergy;
    }

    double energy = 0.0;
    for (int i = 0; i < nthr; ++i) {
        energy += energyArray[i];
    }

    // update number generators
    rng_space.discard(largeNumber * count);
    rng_momentum.discard(largeNumber * count);

    return {energy, accepted_count};
}

double distribute_initial_particles(ParticlesArray& species,
                                    const std::vector<std::unique_ptr<IDistribution>>& distributions,
                                    const Domain& domain, bool co_locate_species) {
    RECORD_TIMER;
    double total_energy = 0.0;

    // Общий seed координат для всех сортов даёт совместную раскладку (co-location).
    // При co_locate_species=false seed зависит от имени сорта, поэтому координаты
    // разных сортов разъезжаются.
    const uint64_t spaceSeed = co_locate_species ? 13 : 13 + hash(species.name(), 1000003);
    LehmerEngine randGenSpace(spaceSeed);
    LehmerEngine randGenPulse(15);

    for (const std::unique_ptr<IDistribution>& dist : distributions) {
        total_energy +=
            add_particles_from_distribution(species, *dist, randGenSpace, randGenPulse, domain, 0.0, true).energy;
    }

    return total_energy;
}

double inject_particles_step(ParticlesArray& species, std::vector<std::unique_ptr<IDistribution>>& distributions,
                             int timestep, const Domain& domain, double dt, bool co_locate_species) {
    RECORD_TIMER;
    if (distributions.empty()) {
        return 0.0;
    }

    double step_energy = 0.0;
    int step_count = 0;

    const int seedOffset = co_locate_species ? 0 : hash(species.name(), 1000003);
    LehmerEngine randGenSpace(13 * 3 * timestep + seedOffset);
    LehmerEngine randGenPulse(hash(species.name(), 20) + 3 * timestep);

    for (const std::unique_ptr<IDistribution>& dist : distributions) {
        // Advance the fractional rate once; count and energy describe the same
        // particles after boundary clipping, including bound-injection drift.
        const InjectedParticles injected =
            add_particles_from_distribution(species, *dist, randGenSpace, randGenPulse, domain, dt, true);
        step_count += injected.count;
        step_energy += injected.energy;
    }

    species.step_state.injection_count += step_count;
    return step_energy;
}
