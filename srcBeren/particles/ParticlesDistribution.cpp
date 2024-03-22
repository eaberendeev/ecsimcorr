#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "service.h"

void ParticlesArray::set_space_distribution() {
    int startX = 0;
    int endX = _world.regionGlob.numCells.x();
    int startY = 0;
    int endY = _world.regionGlob.numCells.y();
    int startZ = 0;
    int endZ = _world.regionGlob.numCells.z();
    // int startY = 0.5*(_world.regionGlob.numCells.y() - int(widthY/Dy));
    // int endY = 0.5*(_world.regionGlob.numCells.y() + int(widthY/Dy));
    // int startZ = 0.5*(_world.regionGlob.numCells.z() - int(widthZ/Dz));
    // int endZ = 0.5*(_world.regionGlob.numCells.z() + int(widthZ/Dz));

    if (initDist == "StrictUniform") {
        if (NumPartPerLine * NumPartPerLine * NumPartPerLine !=
            NumPartPerCell) {
            std::cout << "Wrong value NumPartPerLine! Please use "
                         "NumPartPerCell = NumPartPerLine**3\n";
            exit(-1);
        }
        set_strict_uniform(int3(startX, startY, startZ),
                           int3(endX, endY, endZ));
    }
    if (initDist == "Uniform") {
        RandomGenerator randomGenerator;
        randomGenerator.SetRandSeed(200);
        set_uniform(int3(startX, startY, startZ), int3(endX, endY, endZ),
                    randomGenerator);
    }

    if (initDist == "StrictUniformCircle") {
        if (NumPartPerLine * NumPartPerLine * NumPartPerLine !=
            NumPartPerCell) {
            std::cout << "Wrong value NumPartPerLine! Please use "
                         "NumPartPerCell = NumPartPerLine**3\n";
            exit(-1);
        }
        set_strict_uniform(int3(startX, startY, startZ),
                           int3(endX, endY, endZ));
        set_uniform_circle(int3(startX, startY, startZ),
                           int3(endX, endY, endZ));
    }
    if (initDist == "UniformCircle") {
        RandomGenerator randomGenerator;
        randomGenerator.SetRandSeed(200);
        set_uniform(int3(startX, startY, startZ), int3(endX, endY, endZ),
                    randomGenerator);
        set_uniform_circle(int3(startX, startY, startZ),
                           int3(endX, endY, endZ));
    }

    std::cout << "Distribution " << name << " : " << initDist << "\n";
    return;
}
void ParticlesArray::set_pulse_distribution(
    ThreadRandomGenerator &randomGenerator) {
    auto sigma = temperature;
    double3 pulse;
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            // double sigmaXY = velocity_from_kev(0.5);
            pulse.x() = randomGenerator.Gauss(temperature / sqrt(_mass));
            pulse.y() = randomGenerator.Gauss(temperature / sqrt(_mass));
            // double sigmaZ = velocity_from_kev(1.0);
            pulse.z() = randomGenerator.Gauss(temperature / sqrt(_mass));
            particle.velocity = pulse;
        }
    }
}

void ParticlesArray::set_strict_uniform(int3 start, int3 end) {
    const double3 delta =
        double3(Dx / NumPartPerLine, Dy / NumPartPerLine, Dz / NumPartPerLine);

    Particle particle;
    double x, y, z;

    x = -0.5 * delta.x();
    for (auto i = 0; i < NumPartPerLine * (end.x() - start.x()); ++i) {
        x += delta.x();
        y = -0.5 * delta.y() + Dy * start.y();
        if (!_world.region.in_region(x + Dx * start.x()))
            continue;

        for (auto j = 0; j < NumPartPerLine * (end.y() - start.y()); ++j) {
            y += delta.y();
            z = -0.5 * delta.z() + Dz * start.z();
            for (auto k = 0; k < NumPartPerLine * (end.z() - start.z()); ++k) {
                z += delta.z();

                particle.coord.x() = x + Dx * start.x() - _world.region.origin;
                particle.coord.y() = y;
                particle.coord.z() = z;
                add_particle(particle);
            }
        }
    }
    update_count_in_cell();
}

void ParticlesArray::set_uniform(int3 start, int3 end,
                                 RandomGenerator &randomGenerator) {
    Particle particle;
    double x, y, z;

    for (auto i = 0; i < NumPartPerCell * (end.y() - start.y()) *
                             (end.x() - start.x()) * (end.z() - start.z());
         ++i) {
        x = Dx * (start.x() + randomGenerator.Uniform01() * (end.x() - start.x()));
        y = Dy * (start.y() + randomGenerator.Uniform01() * (end.y() - start.y()));
        z = Dz * (start.z() + randomGenerator.Uniform01() * (end.z() - start.z()));
        if (!_world.region.in_region(x))
            continue;

        particle.coord.x() = x - _world.region.origin;
        particle.coord.y() = y;
        particle.coord.z() = z;
        add_particle(particle);
    }
    update_count_in_cell();
}

void ParticlesArray::set_uniform_circle(int3 start, int3 end) {
    double y, z;
    double center_y = 0.5 * Dy * (end.y() + start.y());
    double center_z = 0.5 * Dz * (end.z() + start.z());
    double radius_y = 0.5 * Dy * (end.y() - start.y());
    double radius_z = 0.5 * Dz * (end.z() - start.z());

    for (int ix = 0; ix < particlesData.size().x(); ++ix) {
        for (int iy = 0; iy < particlesData.size().y(); ++iy) {
            for (int iz = 0; iz < particlesData.size().z(); ++iz) {
                int ip = 0;
                while (ip < countInCell(ix, iy, iz)) {
                    Particle particle = particlesData(ix, iy, iz)[ip];

                    y = particle.coord.y() - center_y;
                    z = particle.coord.z() - center_z;
                    bool outOfCircle = y * y / (radius_y * radius_y) +
                                           z * z / (radius_z * radius_z) >
                                       1.;
                    if (outOfCircle) {
                        delete_particle_runtime(ix, iy, iz, ip);
                    } else {
                        ip++;
                    }
                }
            }
        }
    }
    update_count_in_cell();
}
