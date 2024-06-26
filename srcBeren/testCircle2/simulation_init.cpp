// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <iostream>

#include "recovery.h"
#include "simulation.h"
#include "util.h"

void Simulation::init(Mesh &mesh) {
    mesh.init(domain, parameters);

    set_fields(mesh);
    set_particles();
    make_folders();
}

void Simulation::set_particles() {
    std::vector<ParametersMap> vecParamsMap =
        load_vector_parameters("./PartParams.cfg", "Particles");

    for (const auto &particlesParameters : vecParamsMap) {
        species.emplace_back(particlesParameters, parameters, domain);
    }
    for (auto &sp : species) {
        if (parameters.get_int("StartFromTime") > 0) {
            read_particles_from_recovery(sp);
            std::cout << "Upload " + sp.name() + " success!\n";
            continue;
        }
        if (parameters.get_int("k_particles_reservation") > 0.) {
            for (auto k = 0; k < sp.size(); ++k) {
                sp.particlesData(k).reserve(
                    parameters.get_int("k_particles_reservation") *
                    parameters.get_int("NumPartPerCell"));
            }
        }
            sp.density_on_grid_update();
            std::cout << sp.particlesData.size() << " "
                      << sp.particlesData.capacity() << "\n";
    }
}

void Simulation::set_fields(Mesh &mesh) {
    mesh.fieldJp.set_zero();
    mesh.fieldJe.set_zero();
    mesh.fieldEn.set_zero();
    mesh.fieldEp.set_zero();

    if (parameters.get_int("StartFromTime") > 0) {
        read_fields_from_recovery(mesh.fieldEn, mesh.fieldB);
    } else {
        mesh.set_uniform_field(mesh.fieldEn, 0.0, 0.0, 0.0);
        mesh.set_uniform_field(mesh.fieldBInit,
                               parameters.get_double("BUniform", 0),
                               parameters.get_double("BUniform", 1),
                               parameters.get_double("BUniform", 2));
        mesh.fieldB = mesh.fieldBInit;
    }
}

void Simulation::set_init_particles_distribution_rectangle(ParticlesArray &sp) {
    ThreadRandomGenerator randGenSpace;
    static ThreadRandomGenerator randGenPulse;

        double3 startCoord(0, 0, 0);
        double3 endCoord;
        endCoord.x() = parameters.get_double("Dx") * parameters.get_int("NumCellsX_glob");
        endCoord.y() = parameters.get_double("Dx") * parameters.get_int("NumCellsY_glob");
        endCoord.z() = parameters.get_double("Dz") * parameters.get_int("NumCellsZ_glob");
        double volume = parameters.get_double("Dx") * parameters.get_double("Dy") * parameters.get_double("Dz");
        int numParticles = endCoord.x() * endCoord.y() * endCoord.z() *
                           parameters.get_int("NumPartPerCell") / volume;
                           for(auto& sp: species)
        double initEnergy = sp.add_uniform_rectangle(
            numParticles, double3(sp.temperature, sp.temperature, sp.temperature),
            startCoord, endCoord, randGenSpace, randGenPulse);

    std::cout << "Distribution " << sp.name() << " is done\n";
}
void Simulation::set_init_particles_distribution_cilinder(ParticlesArray &sp) {
    ThreadRandomGenerator randGenSpace;
    static ThreadRandomGenerator randGenPulse;
    double3 center;
    center.x() =
        parameters.get_double("Dx") * parameters.get_int("NumCellsX_glob") / 2;
    center.y() =
        parameters.get_double("Dx") * parameters.get_int("NumCellsY_glob") / 2;
    center.z() =
        parameters.get_double("Dz") * parameters.get_int("NumCellsZ_glob") / 2;
    double rz =
        parameters.get_double("Dz") * parameters.get_int("NumCellsZ_glob") / 2;
    double rr = 30 * parameters.get_double("Dx");
    int numParticles =
        M_PI * rr * rr * (2 * rz) * parameters.get_int("NumPartPerCell");
    sp.injectionEnergy = sp.add_uniform_cilinderZ(
        numParticles, double3(sp.temperature, sp.temperature, sp.temperature), center,
        rr, rz, randGenSpace, randGenPulse);

    std::cout << "Distribution " << sp.name() << " is done\n";
}
