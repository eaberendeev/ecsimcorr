#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "particles_distribution_collection.h"
#include "service.h"
#include "util.h"

void ParticlesArray::initializeDistributions(const nlohmann::json& config) {
    // NB: убедитесь, что xCellSize доступна (член класса или глобальная
    // переменная)
    if (config.contains("distribution") && config["distribution"].is_array()) {
        for (const auto& dist_config : config["distribution"]) {
            if (!dist_config.contains("type")){
                throw std::runtime_error("Distribution type not specified");
            }
            std::string dist_type = dist_config.value("type", "");

            if ((dist_type == "initial" || dist_type == "injection" ||
                 dist_type == "injection_bound") &&
                dist_config.contains("dist_space") &&
                dist_config.contains("dist_pulse")) {
                Distribution dst;
                dst.position = PositionDistributionFactory::create(
                    dist_config["dist_space"]);
                dst.velocity = VelocityDistributionFactory::create(
                    dist_config["dist_pulse"], _mass);
                double dens = 1;
                if (dist_config.contains("density")){
                    dens = dist_config.value("density", 1.0);
                } else{
                    std::cout << "Density not specified, setting to 1.0" << std::endl;
                }
                // Преобразуем double -> int аккуратно и не даём отрицательных
                // значений xCellSize должен быть определён (проверьте)
                double cell_vol =
                    xCellSize * xCellSize * xCellSize;   // <-- must exist
                double calc = static_cast<double>(NumPartPerCell) *
                              dst.position->get_volume() * dens / cell_vol;
                dst.count = std::max(0, static_cast<int>(std::lround(calc)));
                dst.type = dist_type;

                if (dist_type == "initial")
                    initialDistributions.push_back(dst);
                else
                    injectionDistributions.push_back(dst);
            }
        }
    }
}

double ParticlesArray::add_particles(
    ThreadRandomGenerator& rng_space, ThreadRandomGenerator& rng_momentum,
    const std::vector<Distribution>& distributions, const Domain& domain,
    double dt) {

    double energy = 0;

    if (distributions.empty())
        return energy;

    for (const auto& dist : distributions) {
        int count = dist.get_count();
        for (int i = 0; i < count; ++i) {
            double3 position = dist.position->sample(rng_space);
            double3 velocity = dist.velocity->sample(rng_momentum);
            if (dist.get_type() == "injection_bound") {
                position += velocity * dt;
            }
            Particle particle(position, velocity);
            if(particle_boundaries(particle, domain)){
                energy += get_energy_particle(velocity, _mass, _mpw);
                add_particle(particle);
            }
        }
    }
    return energy;
}

double ParticlesArray::distribute_particles(
    const std::vector<Distribution>& distributions, const Domain& domain,
    double timestep, double dt) {
    ThreadRandomGenerator randGenSpace;
    ThreadRandomGenerator randGenPulse;
    randGenSpace.SetRandSeed(13 + 3 * timestep);
    randGenPulse.SetRandSeed(hash(name(), 20) + 3 * timestep);

    return add_particles(randGenSpace, randGenPulse, distributions, domain, dt);
}
