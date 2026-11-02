#include "ParticlesArray.h"

#include <iterator>
#include <vector>

#include "config/speciesconfig.h"
#include "containers.h"
#include "grid/domain.h"
#include "timer.h"

ParticlesArray::ParticlesArray(const nlohmann::json& config, const Domain& domain)
    : parameters_(species_parameters_from_config(config)), storage_(domain), currentOnGrid(domain.size(), 3) {
    diag.clear();
    diag.is_cylinder_geometry = domain.geom.use_cylinder;
    initialize_distributions(config);
}

ParticlesArray* find_species(Species& species, const std::string& name) {
    auto it = species.find(name);
    return (it == species.end()) ? nullptr : it->second.get();
}

const ParticlesArray* find_species(const Species& species, const std::string& name) {
    auto it = species.find(name);
    return (it == species.end()) ? nullptr : it->second.get();
}

void ParticlesArray::prepare() {
    // Step bookkeeping is independent of whether/when diagnostics are read.
    kineticEnergy = get_kinetic_energy();
    diag.clear();
    currentOnGrid.setZero();
}

void ParticlesArray::move(double dt) {
    timer::commonTimer timer(std::source_location::current().function_name(), -1, timer::MeasureUnit::byte);
    int64_t totalParticles = 0;
#pragma omp parallel for schedule(dynamic, 32) reduction(+ : totalParticles)
    for (auto k = 0; k < size(); ++k) {
        totalParticles += std::ssize(cells()(k));
        for (auto& particle : cells()(k)) {
            particle.move(dt);
        }
    }
    timer.flat.m = sizeof(Particle) * totalParticles;
}
