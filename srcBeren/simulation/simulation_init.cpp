// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <iostream>

#include "recovery.h"
#include "simulation.h"
#include "util.h"

void Simulation::init(Mesh &mesh, World &world) {
    mesh.init(world, domain, parameters);

    set_fields(mesh);
    set_particles(world);
    make_folders();
}

// void Simulation::init() {
//     const int timestep = 0;

//     set_fields();

//     mesh.stencil_matI();
//     mesh.stencil_curlE();
//     mesh.stencil_curlB();
//     mesh.stencil_divE();
//     mesh.stencil_matM(parameters.dt());

//     set_particles();

//     make_folders();
//     output_all(timestep);
//     std::cout << "Initialization step has been finished!\n";
// }

void Simulation::set_particles(World &world) {
    std::vector<std::vector<std::string> > stringParams;
    read_params_to_string("Particles", "./PartParams.cfg", stringParams);
    for (const auto &params : stringParams) {
        species.emplace_back(params, world, domain);
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

    if (RECOVERY > 0) {
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
        endCoord.x() = Dx * NumCellsX_glob;
        endCoord.y() = Dx * NumCellsY_glob;
        endCoord.z() = Dz * NumCellsZ_glob;
        double volume = Dx * Dy * Dz;
        int numParticles = endCoord.x() * endCoord.y() * endCoord.z() *
                           NumPartPerCell / volume;
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
    center.x() = Dx * NumCellsX_glob / 2;
    center.y() = Dx * NumCellsY_glob / 2;
    center.z() = Dz * NumCellsZ_glob / 2;
    double rz = Dz * NumCellsZ_glob / 2;
    double rr = 30 * Dx;
    int numParticles = PI * rr * rr * (2 * rz) * NumPartPerCell;
    sp.injectionEnergy = sp.add_uniform_cilinderZ(
        numParticles, double3(sp.temperature, sp.temperature, sp.temperature), center,
        rr, rz, randGenSpace, randGenPulse);

    std::cout << "Distribution " << sp.name() << " is done\n";
}
