#include "particles_distribution_collection.h"

void distribute_uniform_cilinderZ(
    std::vector<Particle>& particles, const int count, const double3& c,
    const double r0, const double z0, ThreadRandomGenerator& randGenSpace) {
    particles.resize(count);
    double rx, ry;
    for (auto& particle : particles) {
        do {
            rx = (1 - 2 * randGenSpace.Uniform01()) * r0;
            ry = (1 - 2 * randGenSpace.Uniform01()) * r0;

        } while (rx * rx + ry * ry > r0 * r0);

        const double x = c.x() + rx;
        const double y = c.y() + ry;
        const double z = c.z() + 0.5*z0 * (1 - 2 * randGenSpace.Uniform01());

        particle.coord = double3(x, y, z);
    }
    return;
}

void distribute_uniform_rectangle(std::vector<Particle>& particles,
                                  const int count, const double3& c,
                                  const double3& length,
                                  ThreadRandomGenerator& randGenSpace) {
    particles.resize(count);
    for (auto& particle : particles) {
        const double x =
            c.x() + 0.5 * length.x() * (1 - 2 * randGenSpace.Uniform01());
        const double y =
            c.y() + 0.5 * length.y() * (1 - 2 * randGenSpace.Uniform01());
        const double z =
            c.z() + 0.5 * length.z() * (1 - 2 * randGenSpace.Uniform01());

        particle.coord = double3(x, y, z);
    }
    return;
}

void distribute_pulse_gauss(std::vector<Particle>& particles,
                            const double3& sigma,
                            ThreadRandomGenerator& randGenPulse) {
    double3 pulse;

    for (auto& particle : particles) {
        do {
            pulse.x() = randGenPulse.Gauss(sigma.x());
            pulse.y() = randGenPulse.Gauss(sigma.y());
            pulse.z() = randGenPulse.Gauss(sigma.z());
        } while (fabs(pulse.x()) > 3 * sigma.x() ||
                 fabs(pulse.y()) > 3 * sigma.y() ||
                 fabs(pulse.z()) > 3 * sigma.z());

        particle.velocity = pulse;
    }
}
