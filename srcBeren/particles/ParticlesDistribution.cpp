#include <omp.h>

#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "particles_distribution_collection.h"
#include "random.h"

void ParticlesArray::initialize_distributions(const nlohmann::json& config) {
    const double cell_volume = domain_.cell_size().elements_product();
    // Используем фабричный метод для создания всех распределений
    std::vector all_distributions =
        DistributionFactory::createFromConfig(config, cell_volume, NumPartPerCell, mass_, mpw_);

    // Разделяем по типам
    for (auto& dist : all_distributions) {
        std::string type = dist->get_type();
        if (type == "initial") {
            initialDistributions_.push_back(std::move(dist));
        } else if (type == "injection" || type == "injection_bound") {
            injectionDistributions_.push_back(std::move(dist));
        }
    }
}

double ParticlesArray::add_particles_from_distribution(IDistribution& dist, LehmerEngine& rng_space,
                                                       LehmerEngine& rng_momentum, const Domain& domain, double dt,
                                                       bool check_boundaries = true) {
    const int count = dist.get_count_to_inject();
    RECORD_TIMER_PARAMS(count);

    const bool useThreads = (count > 15'000);

    std::vector<std::atomic<int64_t>> locks(useThreads ? domain.total_size() : 0);
    if (useThreads) {
#pragma omp parallel for
        for (int i = 0; i < std::ssize(locks); ++i) {
            locks[i].store(0);
        }
    }

    // run-to run reproducible energy result
    constexpr int64_t largeNumber = 32;
    const int nthr = useThreads ? omp_get_max_threads() : 1;
    double energyArray[nthr];
    std::fill_n(energyArray, nthr, 0.0);
#pragma omp parallel if (useThreads)
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
                threadLocalEnergy += dist.get_energy(velocity);

                if (useThreads) {
                    const Vector3I cell_id = domain_.get_cell_index(particle.coord);
                    std::atomic<int64_t>& lock = locks.at(domain_.sind(cell_id.x(), cell_id.y(), cell_id.z()));
                    int64_t expected = 0;
                    while (!lock.compare_exchange_strong(expected, 1)) {
                        expected = 0;
                        asm volatile("nop;nop;nop;nop;nop;nop;nop;");
                    }
                    particlesData(cell_id.x(), cell_id.y(), cell_id.z()).push_back(particle);
                    lock.store(0);
                } else {
                    add_particle(particle);
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

    return energy;
}

double ParticlesArray::distribute_initial_particles(const std::vector<std::unique_ptr<IDistribution>>& distributions,
                                                    const Domain& domain) {
    RECORD_TIMER;
    double total_energy = 0.0;

    LehmerEngine randGenSpace(13);
    LehmerEngine randGenPulse(15);

    for (auto& dist : distributions) {
        total_energy += add_particles_from_distribution(*dist, randGenSpace, randGenPulse, domain, 0.0, true);
    }

    return total_energy;
}

double ParticlesArray::inject_particles_step(std::vector<std::unique_ptr<IDistribution>>& distributions, int timestep,
                                             const Domain& domain, double dt) {
    RECORD_TIMER;
    if (distributions.empty()) {
        return 0.0;
    }

    double step_energy = 0.0;
    int step_count = 0;

    LehmerEngine randGenSpace(13 * 3 * timestep);
    LehmerEngine randGenPulse(hash(name(), 20) + 3 * timestep);

    for (auto& dist : distributions) {
        step_count += dist->get_count_to_inject();
        step_energy += add_particles_from_distribution(*dist, randGenSpace, randGenPulse, domain, dt, true);
    }

    diag.injection_count += step_count;
    return step_energy;
}

// Метод для добавления распределения во время выполнения
void ParticlesArray::add_distribution(const nlohmann::json& config, const std::string& type) {
    const double cell_volume = domain_.cell_size().elements_product();

    auto dist = DistributionFactory::create(config, type, cell_volume, NumPartPerCell, mass_, mpw_);

    if (type == "initial") {
        initialDistributions_.push_back(std::move(dist));
    } else if (type == "injection" || type == "injection_bound") {
        injectionDistributions_.push_back(std::move(dist));
    }
}
