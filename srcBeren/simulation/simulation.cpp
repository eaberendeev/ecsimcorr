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
#include "Vec.h"
#include "World.h"
#include "collision.h"

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

    // set_fields(domain,parameters);
    // set_species(domain,parameters);



            Region regionGlob;

    Region region = regionGlob;
    /// Make modeling area
    World world(regionGlob, region);

    ///	Make mesh
    Mesh mesh(world);
    BinaryCollider collider(n0,15);
    std::vector<std::vector<std::string> > stringParams;
    ///// Make particles
    std::vector<ParticlesArray> species;
    read_params_to_string("Particles", "./PartParams.cfg", stringParams);
    for (const auto &params : stringParams) {
        species.emplace_back(params, world);
    }

    Writer writer(world, mesh, species);

    writer.output(0.0, StartTimeStep);

    std::map<std::string, FILE *> fDiagParticles;

    for (auto &sp : species) {
        fDiagParticles[sp.name] = fopen((sp.name + ".dat").c_str(), "w");
    }

    Timer globalTimer("globalFunctions.time");
    double energyP, energyPn;
    for (auto timestep = StartTimeStep + 1; timestep <= 2000000; ++timestep) {
        globalTimer.start("Total");

        mesh.prepare();
        for (auto &sp : species) {
            sp.delete_bounds();
            sp.inject_particles(timestep);
            sp.prepare(timestep);   // save start coord for esirkepov current
        }

        // collider.collide_particles(species, NumPartPerCell, Dt);

        globalTimer.start("densityCalc");
        for (auto &sp : species) {
            sp.density_on_grid_update();
        }
        globalTimer.finish("densityCalc");

        collect_charge_density(mesh.chargeDensityOld, species);

        globalTimer.start("particles1");
        for (auto &sp : species) {
            //sp.move(1. * Dt);   // +++ x_{n+1/2} -> x_{n+1}
            sp.move_and_calc_current(0.5 * Dt);   //  +++ x_n -> x_{n+1/2}
            sp.update_cells();
            // +++ get J(x_{n+1/2},v_n)_predict
            sp.predict_current(mesh.fieldB, mesh.fieldJp);
            // +++ get Lgg'(x_{n+1/2})
            sp.get_L(mesh);
        }
        globalTimer.finish("particles1");

        mesh.glue_Lmat_bound();

        globalTimer.start("stencilLmat");
        mesh.stencil_Lmat();
        globalTimer.finish("stencilLmat");

        mesh.make_periodic_border_with_add(mesh.fieldJp);

        globalTimer.start("FieldsPredict");
        // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})). mesh consist En,
        // En+1_predict
        mesh.predictE();
        globalTimer.finish("FieldsPredict");

        globalTimer.start("particles2");
        energyP = energyPn = 0.;
        for (auto &sp : species) {
            energyP += sp.get_kinetic_energy();
            // +++ get v'_{n+1} from v_{n} and E'_{n+1}
            sp.predict_velocity(mesh);
            // +++ x_{n+1/2} -> x_{n+1}
            //sp.move(0.5 * Dt);
            sp.move_and_calc_current(0.5 * Dt);

            sp.update_cells();
            mesh.make_periodic_border_with_add(sp.currentOnGrid);
            sp.currentOnGrid.data() *= 0.5;
        }
        globalTimer.finish("particles2");

        collect_current(mesh.fieldJe, species);

        globalTimer.start("FieldsCorr");
        // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final
        mesh.correctE();
        globalTimer.finish("FieldsCorr");

        globalTimer.start("particles3");
        for (auto &sp : species) {
            sp.correctv(mesh);
            energyPn += sp.get_kinetic_energy();   // / sp.mpw(0);
        }
        for (auto &sp : species) {
            sp.density_on_grid_update();
        }
        globalTimer.finish("particles3");

        globalTimer.start("computeB");
        mesh.computeB();
        globalTimer.finish("computeB");

        collect_charge_density(mesh.chargeDensity, species);
        std::cout << mesh.chargeDensity.data().norm()
                  << " norm mesh.chargeDensity \n";

        auto divJ = mesh.divE * mesh.fieldJe.data();

        auto delta =
            (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (Dt) +
            divJ;
        std::cout << delta.norm() << " norm drho / Dt - divJ \n";
        globalTimer.start("Output");
        writer.output((energyPn - energyP), timestep);
        globalTimer.finish("Output");
        //// DAMPING //////

        mesh.fieldB.data() -= mesh.fieldBInit.data();
        //damping_fields(mesh.fieldEn, mesh.fieldB, region);
        damping_fields_circleXY(mesh.fieldEn, mesh.fieldB, region);
        mesh.fieldB.data() += mesh.fieldBInit.data();

        /// END DAMPING ///////
        globalTimer.finish("Total");
        globalTimer.write(timestep);
        for (auto &sp : species) {
            sp.write_particles_to_recovery(timestep);
        }
        mesh.write_fields_to_recovery(timestep);
    }
}

void Simulation::collect_charge_density(Field3d &field,
                            const std::vector<ParticlesArray> &species) {
    field.clear();
    for (const auto &sp : species) {
#pragma omp parallel for collapse(3)
        for (auto i = 0; i < sp.densityOnGrid.size().x(); i++) {
            for (auto j = 0; j < sp.densityOnGrid.size().y(); j++) {
                for (auto k = 0; k < sp.densityOnGrid.size().z(); k++) {
                    field(i, j, k, 0) += sp.densityOnGrid(i, j, k);
                }
            }
        }
    }
}

void Simulation::collect_current(Field3d &field,
                                 const std::vector<ParticlesArray> &species) {
    field.clear();
    for (const auto &sp : species) {
        field.data() += sp.currentOnGrid.data();
    }
}
