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
    std::vector<std::vector<std::string> > stringParams;
    read_params_to_string("Particles", "./PartParams.cfg", stringParams);
    for (const auto &particlesParameters : stringParams) {
        species.emplace_back(particlesParameters, parameters, domain);
    }
    for (auto &sp : species) {
        if (RECOVERY > 0) {
            read_particles_from_recovery(sp);
            std::cout << "Upload " + sp.name() + " success!\n";
            continue;
        }
        if (k_particles_reservation > 0.) {
            for (auto k = 0; k < sp.size(); ++k) {
                sp.particlesData(k).reserve(
                    int(k_particles_reservation * NumPartPerCell));
            }
        }
        set_init_particles_distribution_rectangle(sp);
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

    mesh.set_uniform_field(mesh.fieldEn, 0.0, 0.0, 0.0);
    mesh.set_uniform_field(mesh.fieldBInit,
                           parameters.get_double("BUniform", 0),
                           parameters.get_double("BUniform", 1),
                           parameters.get_double("BUniform", 2));
    mesh.fieldB = mesh.fieldBInit;
}

void Simulation::set_init_particles_distribution_rectangle(ParticlesArray &sp) {
    ThreadRandomGenerator randGenSpace;
    static ThreadRandomGenerator randGenPulse;

    double3 startCoord(0, 0, 0);
    double3 endCoord;
    endCoord.x() = Dx * NumCellsX_glob;
    endCoord.y() = Dx * NumCellsY_glob;
    endCoord.z() = Dz * NumCellsZ_glob;
    double volume = Dx * Dy * Dz;
    double Txy = sp.temperature;   // sqrt(0.01 / 512.);
    double Tz = sp.temperature;    // sqrt(0.05 / 512.);
    int numParticles =
        endCoord.x() * endCoord.y() * endCoord.z() * NumPartPerCell / volume;
    //for (auto &sp : species) {
        double initEnergy = sp.add_uniform_rectangle(
            numParticles,
            double3(Txy, Txy, Tz), startCoord,
            endCoord, randGenSpace, randGenPulse);
        std::cout << "Distribution " << sp.get_total_num_of_particles() << " "
                  << sp.name() << " is done\n";
    //}
}
