#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "particles_distribution_collection.h"
#include "service.h"
#include "util.h"


void ParticlesArray::distribute_particles(const ParametersMap& parameters,
                                     double timestep) {
    ThreadRandomGenerator randGenSpace;
    ThreadRandomGenerator randGenPulse;
    randGenSpace.SetRandSeed(13 + 3 * timestep);
    randGenSpace.SetRandSeed(hash(name(),20) + 3 * timestep);

    std::vector<Particle> particles =
        distribute_particles_in_space(parameters, randGenSpace);
    injectionEnergy =
        distribute_particles_pulse(particles, parameters, randGenPulse);
    add_particles(particles);
}

std::vector<Particle> ParticlesArray::distribute_particles_in_space(
    const ParametersMap& parameters, ThreadRandomGenerator& randGenSpace) {
    std::vector<Particle> particles;

    if (distSpace[0] == "UniformCylZ_cx_cy_cz_rr_rz") {
        double3 center;
        center.x() = stod(distSpace[1]);
        center.y() = stod(distSpace[2]);
        center.z() = stod(distSpace[3]);
        double rr = stod(distSpace[4]);
        double rz = stod(distSpace[5]);

        const double cellVolume = xCellSize * yCellSize * zCellSize;
        const int NumPartPerCell = parameters.get_int("NumPartPerCell");

        int count = M_PI * rr * rr * (2 * rz) * NumPartPerCell / cellVolume;
        std::cout << distType << std::endl;
        if (distType == "INJECTION") {
            const double dt = parameters.get_double("Dt");
            count = count * dt / parameters.get_double("Tau");
            std::cout << "Particles per step: " << count << std::endl;
        }

        distribute_uniform_cilinderZ(particles, count, center, rr, rz,
                                     randGenSpace);
    } else if (distSpace[0] == "Uniform_cx_cy_cz_lx_ly_lz") {
        double3 center, length;
        center.x() = stod(distSpace[1]);
        center.y() = stod(distSpace[2]);
        center.z() = stod(distSpace[3]);
        length.x() = stod(distSpace[4]);
        length.y() = stod(distSpace[5]);
        length.z() = stod(distSpace[6]);

        const double cellVolume = xCellSize * yCellSize * zCellSize;
        const int NumPartPerCell = parameters.get_int("NumPartPerCell");

        int count = length.x() * length.y() * length.z() * NumPartPerCell / cellVolume;
        std::cout << distType << std::endl;
        if (distType == "INJECTION") {
            const double dt = parameters.get_double("Dt");
            count = count * dt / parameters.get_double("Tau");
            std::cout << "Particles per step: " << count << std::endl;
        }

        distribute_uniform_rectangle(particles, count, center, length,
                                     randGenSpace);
    } else if (distSpace[0] == "None") {
        std::cout << "Particles was not added" << std::endl;
    } else {
        std::cout << "Error: unknown distribution space type" << std::endl;
        exit(1);
    }
    return particles;
}

double ParticlesArray::distribute_particles_pulse(
    std::vector<Particle>& particles, const ParametersMap& parameters,
    ThreadRandomGenerator& randGenPulse) {
    if (distPulse[0] == "Gauss") {
        double3 sigma = (1.0 / sqrt(_mass))* temperature ;
        distribute_pulse_gauss(particles, sigma, randGenPulse);
    } else if (distPulse[0] == "None") {
    } else {
        std::cout << "Error: unknown distribution pulse type" << std::endl;
        exit(1);
    }

    return get_energy_particles(particles, _mass, _mpw);
}
