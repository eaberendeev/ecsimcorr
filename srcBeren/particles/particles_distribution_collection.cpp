#include "particles_distribution_collection.h"

void distribute_uniform_cilinderZ(
    std::vector<Particle>& particles, const int count, const double3& c,
    const double r0, const double rz, ThreadRandomGenerator& randGenSpace) {
    particles.resize(count);
    double rx, ry;
    for (auto& particle : particles) {
        do {
            rx = (1 - 2 * randGenSpace.Uniform01()) * r0;
            ry = (1 - 2 * randGenSpace.Uniform01()) * r0;

        } while (rx * rx + ry * ry > r0 * r0);

        const double x = c.x() + rx;
        const double y = c.y() + ry;
        const double z = c.z() + rz * (1 - 2 * randGenSpace.Uniform01());

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
            c.x() + length.x() * (1 - 2 * randGenSpace.Uniform01());
        const double y =
            c.y() + length.y() * (1 - 2 * randGenSpace.Uniform01());
        const double z =
            c.z() + length.z() * (1 - 2 * randGenSpace.Uniform01());

        particle.coord = double3(x, y, z);
    }
    return;
}

#include <algorithm>
#include <stdexcept>

const double PULSE_LIMIT_FACTOR = 3.0;

void generate_gaussian_pulse(double3& pulse, const double3& sigma,
                             ThreadRandomGenerator& randGenPulse) {
    do {
        pulse.x() = randGenPulse.Gauss(sigma.x());
        pulse.y() = randGenPulse.Gauss(sigma.y());
        pulse.z() = randGenPulse.Gauss(sigma.z());
    } while (fabs(pulse.x()) > PULSE_LIMIT_FACTOR * sigma.x() ||
             fabs(pulse.y()) > PULSE_LIMIT_FACTOR * sigma.y() ||
             fabs(pulse.z()) > PULSE_LIMIT_FACTOR * sigma.z());
}

void distribute_pulse_gauss(std::vector<Particle>& particles,
                            const double3& sigma,
                            ThreadRandomGenerator& randGenPulse) {
    std::for_each(particles.begin(), particles.end(), [&](Particle& particle) {
        double3 pulse;
        generate_gaussian_pulse(pulse, sigma, randGenPulse);
        particle.velocity =
            pulse;   // Only modify velocity, keep other fields intact
    });
}


void distribute_pulse_sin(std::vector<Particle>& particles, const double vx,
                          const double period, const double3& sigma,
                          ThreadRandomGenerator& randGenPulse) {
    if (period <= 0) {
        std::cout<< "Length must be positive\n";
        exit(1);
    }

    // std::for_each(particles.begin(), particles.end(), [&](Particle& particle) {
    //     double3 pulse;
    //     generate_gaussian_pulse(pulse, sigma, randGenPulse);
    //     particle.velocity = pulse;
    //     particle.velocity.x() +=
    //         vx *
    //         sin(2 * M_PI *
    //             (particle.coord.x() / period + particle.coord.y() / period));
    //     particle.velocity.y() +=
    //         vx *
    //         sin(2 * M_PI *
    //             (particle.coord.x() / period + particle.coord.y() / period));
    // });

    std::for_each(particles.begin(), particles.end(), [&](Particle& particle) {
        double3 pulse;
        generate_gaussian_pulse(pulse, sigma, randGenPulse);
        particle.velocity = pulse;
        particle.velocity.y() +=
            vx *
            sin(2 * M_PI *
                (particle.coord.y() / period));
    });
}
