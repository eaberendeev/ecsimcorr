#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "util.h"
#include "service.h"

double ParticlesArray::add_uniform_cilinderZ(
    const int numParts, const double3& temperature, const double3& c,
    const double r0, const double z0, ThreadRandomGenerator& randGenSpace,
    ThreadRandomGenerator& randGenPulse) {
    Particle particle;
    double3 coord, pulse;
    double rx, ry, x, y, z;
    double energy = 0;
    for (int k = 0; k < numParts; k++) {
        do {
            rx = (1 - 2 * randGenSpace.Uniform01()) * r0;
            ry = (1 - 2 * randGenSpace.Uniform01()) * r0;

        } while (rx * rx + ry * ry > r0 * r0);

        x = c.x() + rx;
        y = c.y() + ry;
        z = c.z() + z0 * (1 - 2 * randGenSpace.Uniform01());
        particle.coord = double3(x, y, z);

        pulse.x() = randGenPulse.Gauss(temperature.x() / sqrt(_mass));
        pulse.y() = randGenPulse.Gauss(temperature.y() / sqrt(_mass));
        pulse.z() = randGenPulse.Gauss(temperature.z() / sqrt(_mass));
        particle.velocity = pulse;
        energy += get_energy_particle(particle.velocity, _mass, _mpw);

        add_particle(particle);
    }
    update_count_in_cell();
    return energy;
}

double ParticlesArray::add_uniform_rectangle(
    const int numParts, const double3& temperature, const double3& r,
    const double3& sizeL, ThreadRandomGenerator& randGenSpace,
    ThreadRandomGenerator& randGenPulse) {
    Particle particle;
    double3 coord, pulse;
    double x, y, z;
    double energy = 0;

    for (int k = 0; k < numParts; k++) {
        x = r.x() + sizeL.x() * randGenSpace.Uniform01();
        y = r.y() + sizeL.y() * randGenSpace.Uniform01();
        z = r.z() + sizeL.z() * randGenSpace.Uniform01();
        particle.coord = double3(x, y, z);
        do {
            pulse.x() = randGenPulse.Gauss(temperature.x() / sqrt(_mass));
            pulse.y() = randGenPulse.Gauss(temperature.y() / sqrt(_mass));
            pulse.z() = randGenPulse.Gauss(temperature.z() / sqrt(_mass));
        } while (fabs(pulse.x()) > 3 * temperature.x() ||
                 fabs(pulse.y()) > 3 * temperature.y() ||
                 fabs(pulse.z()) > 3 * temperature.z());

        particle.velocity = pulse;
        energy += get_energy_particle(particle.velocity, _mass, _mpw);

        add_particle(particle);
    }
    update_count_in_cell();
    return energy;
}
