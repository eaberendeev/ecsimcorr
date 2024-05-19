#include "ParticlesArray.h"

void ParticlesArray::inject_particles(const int timestep) {
    static ThreadRandomGenerator randGenSpace;
    static ThreadRandomGenerator randGenPulse;

    double3 center;
    center.x() = Dx * NumCellsX_glob / 2;
    center.y() = Dx * NumCellsY_glob / 2;
    center.z() = Dz * NumCellsZ_glob / 2;
    double rz = Dz * NumCellsZ_glob / 2;
    double rr = 30 * Dx;
    int ParticlesPerStep =
        Dt * PI * rr * rr * (2 * rz) * NumPartPerCell / (Tau * Dx * Dy * Dz);
    randGenSpace.SetRandSeed(100 + timestep);
    injectionEnergy = add_uniform_cilinder(ParticlesPerStep, rr, rz, center,
                                           randGenSpace, randGenPulse);
}

double ParticlesArray::add_uniform_cilinder(
    int numParts, double r0, double z0, double3 c,
    ThreadRandomGenerator& randGenSpace, ThreadRandomGenerator& randGenPulse) {
    Particle particle;
    auto sigma = temperature;
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

        pulse.x() = randGenPulse.Gauss(sigma / sqrt(_mass));
        pulse.y() = randGenPulse.Gauss(sigma / sqrt(_mass));
        pulse.z() = randGenPulse.Gauss(sigma / sqrt(_mass));
        particle.velocity = pulse;
        energy += get_energy_particle(particle.velocity, _mass, _mpw);

        add_particle(particle);
    }
    update_count_in_cell();
    return energy;
}

double ParticlesArray::add_uniform_line(int numParts, double3 r, double3 sizeL,
                                      ThreadRandomGenerator& randGenSpace,
                                      ThreadRandomGenerator& randGenPulse) {
    Particle particle;
    auto sigma = temperature;
    double3 coord, pulse;
    double x, y, z;
    double energy = 0;

    for (int k = 0; k < numParts; k++) {
        x = r.x() + sizeL.x() * randGenSpace.Uniform01();
        y = r.y() + sizeL.y() * randGenSpace.Uniform01();
        z = r.z() + sizeL.z() * randGenSpace.Uniform01();
        particle.coord = double3(x, y, z);
        do {
            pulse.x() = randGenPulse.Gauss(sigma / sqrt(_mass));
            pulse.y() = randGenPulse.Gauss(sigma / sqrt(_mass));
            pulse.z() = randGenPulse.Gauss(sigma / sqrt(_mass));
        } while (fabs(pulse.x()) > 3 * sigma || fabs(pulse.y()) > 3 * sigma ||
                 fabs(pulse.z()) > 3 * sigma);

        particle.velocity = pulse;
        energy += get_energy_particle(particle.velocity, _mass, _mpw);

        add_particle(particle);
    }
    update_count_in_cell();
    return energy;
}
