// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file
#include <fstream>
#include <iostream>
#include "containers.h"
#include "ParticlesArray.h"
#include "output_util.h"

// FIELDS

void read_fields_from_recovery(Field3d& fieldE, Field3d& fieldB) {
    std::string fname = "..//Recovery//Fields//FieldE.backup";
    read_field_from_file(fname, fieldE);
    fname = "..//Recovery//Fields//FieldB.backup";
    read_field_from_file(fname, fieldB);
}

void write_fields_to_recovery(const Field3d& fieldE, const Field3d& fieldB,
                              const int timestep, const int recoveryInterval) {
    if (recoveryInterval < 0) {
        return;
    }
    if (timestep % recoveryInterval != 0) {
        return;
    }

    std::cout << "Backup mesh in " << timestep << "\n";
    std::string fname = ".//Recovery//Fields//FieldE.backup";
    write_field_to_file(fname, fieldE);
    fname = ".//Recovery//Fields//FieldB.backup";
    write_field_to_file(fname, fieldB);
}

// PARTICLES
void write_particles_to_recovery(const ParticlesArray& particles,
                                 const int timestep,
                                 const int recoveryInterval) {
    if (recoveryInterval < 0) {
        return;
    }
    if (timestep % recoveryInterval != 0) {
        return;
    }

    std::cout << "Backup " + particles.name() + " in " << timestep << "\n";
    Particle particle;
    std::ofstream file_bin(".//Recovery//Particles//" + particles.name() + "//" +
                               particles.name() + ".backup",
                           std::ios::out | std::ios::binary);
    int n_particles = 0;
    for (auto cell = 0; cell < particles.size(); ++cell) {
        n_particles += particles.particlesData(cell).size();
    }

    file_bin.write((char*) &n_particles, sizeof(n_particles));
    for (auto cell = 0; cell < particles.size(); ++cell) {
        int cellSize = particles.particlesData(cell).size();
        if (cellSize > 0) {
            file_bin.write((char*) &particles.particlesData(cell)[0],
                           cellSize * sizeof(particles.particlesData(cell)[0]));
        }
    }
    file_bin.close();
}

void read_particles_from_recovery(
    ParticlesArray& particles) {
    Particle particle;
    std::ifstream file_bin("..//Recovery//Particles//" + particles.name() + "//" +
                               particles.name() + ".backup",
                           std::ios::in | std::ios::binary);
    int n_particles;
    file_bin.read((char*) &n_particles, sizeof(n_particles));
    for (int i = 0; i < n_particles; i++) {
        file_bin.read((char*) &particle, sizeof(Particle));
        particles.add_particle(particle);
    }
    file_bin.close();
}
