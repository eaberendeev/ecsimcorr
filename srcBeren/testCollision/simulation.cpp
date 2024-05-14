// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation.h"

#include <iostream>
#include <map>
#include <string>

#include "Damping.h"
#include "Diagnostic.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "Read.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "recovery.h"

Simulation::Simulation(int argc, char **argv)
    : parameters(load_parameters("./SysParams.cfg"))
{
    static char help[] = "Plasma simulation.\n\n";
    std::cout << help;
    parameters.print();


}

void Simulation::make_all() {
    bounds.setBounds(Bounds::BoundValues(BoundType::PERIODIC, BoundType::PERIODIC,
                               BoundType::PERIODIC),
                     Bounds::BoundValues(BoundType::PERIODIC, BoundType::PERIODIC,
                               BoundType::PERIODIC));
    domain.setDomain(parameters, bounds);
    const double dt  = parameters.get_double("Dt");

//    BinaryCollider collider(n0,15);

    init(mesh);

    Timer globalTimer("globalFunctions.time");
    int timestep = 0;
    output_all(mesh, species, timestep);
    static std::mt19937 generator;
    BinaryCollider collider;

    for (auto timestep = StartTimeStep + 1; timestep <= parameters.get_int("LastTimestep"); ++timestep) {
        globalTimer.start("Total");

        mesh.prepare();

        for (auto &sp : species) {
            sp.prepare(timestep);   // save start coord for esirkepov current
        }
        globalTimer.start("densityCalc");
        for (auto &sp : species) {
            sp.density_on_grid_update();
        }
        globalTimer.finish("densityCalc");

        //collider.collide_particles(species, NumPartPerCell, Dt);
        //for (auto &sp : species) {
            collider.collide_same_sort_binary(species, Dt);
            collider.collide_ion_electron_binary(species, Dt);
            //}
            globalTimer.start("particles1");
            for (auto &sp : species) {
                sp.move(1. * Dt);   // +++ x_{n+1/2} -> x_{n+1}
                sp.update_cells(domain);
                std::cout << "Moved " << sp.get_total_num_of_particles() << " "
                          << sp.name() << "\n";
        }
        globalTimer.finish("particles1");

        /// END DAMPING ///////
        globalTimer.finish("Total");
        globalTimer.write(timestep, 1);
        output_all(mesh, species, timestep);
    }
}
