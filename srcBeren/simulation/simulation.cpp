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

void Simulation::inject_particles(const int timestep) {
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
    for(auto& sp : species){
        randGenSpace.SetRandSeed(100 + timestep);
        sp.injectionEnergy = sp.add_uniform_cilinderZ(
            ParticlesPerStep, double3(sp.temperature, sp.temperature, sp.temperature),
            center, rr, rz, randGenSpace, randGenPulse);
    }
}

void Simulation::make_all() {
    bounds.setBounds(Bounds::BoundValues(BoundType::PERIODIC, BoundType::PERIODIC,
                               BoundType::PERIODIC),
                     Bounds::BoundValues(BoundType::PERIODIC, BoundType::PERIODIC,
                               BoundType::PERIODIC));
    domain.setDomain(parameters, bounds);
    const double dt  = parameters.get_double("Dt");
    Region regionGlob(parameters);
    Region region = regionGlob;
    /// Make modeling area
    World world(regionGlob, region);

    BinaryCollider collider(n0,15);

    init(mesh, world);
    Writer writer(world, mesh, species, domain, parameters);

    writer.output(0.0, parameters, StartTimeStep);

    // std::map<std::string, FILE *> fDiagParticles;

    // for (auto &sp : species) {
    //     fDiagParticles[sp.name] = fopen((sp.name + ".dat").c_str(), "w");
    // }

    Timer globalTimer("globalFunctions.time");
    double energyP, energyPn;
    for (auto timestep = StartTimeStep + 1; timestep <= 2000000; ++timestep) {
        globalTimer.start("Total");

        mesh.prepare();
        inject_particles(timestep);

        for (auto &sp : species) {
            sp.delete_bounds();
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
            sp.update_cells(domain);
            // +++ get J(x_{n+1/2},v_n)_predict
            sp.predict_current(mesh.fieldB, mesh.fieldJp, domain);
            // +++ get Lgg'(x_{n+1/2})
            sp.get_L(mesh, domain, dt);
        }
        globalTimer.finish("particles1");

        mesh.glue_Lmat_bound();

        globalTimer.start("stencilLmat");
        mesh.stencil_Lmat(domain);
        globalTimer.finish("stencilLmat");

        mesh.make_periodic_border_with_add(mesh.fieldJp);

        globalTimer.start("FieldsPredict");
        // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})). mesh consist En,
        // En+1_predict
        mesh.predictE(dt);
        globalTimer.finish("FieldsPredict");

        globalTimer.start("particles2");
        energyP = energyPn = 0.;
        for (auto &sp : species) {
            energyP += sp.get_kinetic_energy();
            // +++ get v'_{n+1} from v_{n} and E'_{n+1}
            sp.predict_velocity(mesh, domain);
            // +++ x_{n+1/2} -> x_{n+1}
            //sp.move(0.5 * Dt);
            sp.move_and_calc_current(0.5 * Dt);

            sp.update_cells(domain);
            mesh.make_periodic_border_with_add(sp.currentOnGrid);
            sp.currentOnGrid.data() *= 0.5;
        }
        globalTimer.finish("particles2");

        collect_current(mesh.fieldJe, species);

        globalTimer.start("FieldsCorr");
        // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final
        mesh.correctE(dt);
        globalTimer.finish("FieldsCorr");

        globalTimer.start("particles3");
        for (auto &sp : species) {
            sp.correctv(mesh, domain);
            energyPn += sp.get_kinetic_energy();   // / sp.mpw(0);
        }
        for (auto &sp : species) {
            sp.density_on_grid_update();
        }
        globalTimer.finish("particles3");

        globalTimer.start("computeB");
        mesh.computeB(mesh.fieldE, mesh.fieldEn, mesh.fieldB, parameters.get_double("Dt"));
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
        writer.output((energyPn - energyP), parameters, timestep);
        globalTimer.finish("Output");
        //// DAMPING //////

        mesh.fieldB.data() -= mesh.fieldBInit.data();
        //damping_fields(mesh.fieldEn, mesh.fieldB, region);
        damping_fields_circleXY(mesh.fieldEn, mesh.fieldB, domain, parameters);
        mesh.fieldB.data() += mesh.fieldBInit.data();

        /// END DAMPING ///////
        globalTimer.finish("Total");
        globalTimer.write(timestep, 1);
        output_all(mesh, species, timestep);
    }
}

void Simulation::collect_charge_density(Field3d &field,
                            const std::vector<ParticlesArray> &species) {
    field.set_zero();
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
    field.set_zero();
    for (const auto &sp : species) {
        field.data() += sp.currentOnGrid.data();
    }
}
