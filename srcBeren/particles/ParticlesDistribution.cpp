#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "particles_distribution_collection.h"

void ParticlesArray::initialize_distributions(const nlohmann::json& config) {
    const double cell_volume = domain_.cell_size().elements_product();
    // Используем фабричный метод для создания всех распределений
    auto all_distributions = DistributionFactory::createFromConfig(config, cell_volume, NumPartPerCell, mass_, mpw_);

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

double ParticlesArray::add_particles_from_distribution(IDistribution& dist, ThreadRandomGenerator& rng_space,
                                                       ThreadRandomGenerator& rng_momentum, const Domain& domain,
                                                       double dt, bool check_boundaries = true) {
    double energy = 0.0;
    int count = dist.get_count_to_inject();

    for (int i = 0; i < count; ++i) {
        Vector3R position = dist.sample_position(rng_space);
        Vector3R velocity = dist.sample_velocity(rng_momentum);

        if (dist.is_bound_injection()) {
            position += velocity * dt;
        }

        Particle particle(position, velocity);

        if (!check_boundaries || domain.contains(position)) {
            energy += dist.get_energy(velocity);
            add_particle(particle);
        }
    }

    return energy;
}

double ParticlesArray::distribute_initial_particles(const std::vector<std::unique_ptr<IDistribution>>& distributions,
                                                    const Domain& domain) {
    double total_energy = 0.0;

    ThreadRandomGenerator randGenSpace;
    ThreadRandomGenerator randGenPulse;
    randGenSpace.SetRandSeed(13);
    randGenPulse.SetRandSeed(15);

    for (auto& dist : distributions) {
        total_energy += add_particles_from_distribution(*dist, randGenSpace, randGenPulse, domain, 0.0, true);
    }

    return total_energy;
}

double ParticlesArray::inject_particles_step(std::vector<std::unique_ptr<IDistribution>>& distributions, int timestep,
                                             const Domain& domain, double dt) {
    if (distributions.empty()) {
        return 0.0;
    }

    double step_energy = 0.0;

    static ThreadRandomGenerator randGenSpace;
    static ThreadRandomGenerator randGenPulse;
    randGenSpace.SetRandSeed(13 + 3 * timestep);
    randGenPulse.SetRandSeed(hash(name(), 20) + 3 * timestep);

    for (auto& dist : distributions) {
        step_energy += add_particles_from_distribution(*dist, randGenSpace, randGenPulse, domain, dt, true);
    }

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
