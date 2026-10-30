#include "ParticlesArray.h"
#include "Shape.h"
#include "containers.h"
#include "grid/domain.h"
#include "memory.h"

void ParticlesArray::move(double dt) {
    timer::commonTimer timer(std::source_location::current().function_name(), -1, timer::MeasureUnit::byte);
    int64_t totalParticles = 0;
#pragma omp parallel for schedule(dynamic, 32) reduction(+ : totalParticles)
    for (auto k = 0; k < size(); ++k) {
        totalParticles += std::ssize(particlesData(k));
        for (auto& particle : particlesData(k)) {
            particle.move(dt);
        }
    }
    timer.flat.m = sizeof(Particle) * totalParticles;
}
