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

// TO DO:: call function set_init_particles_distribution_cilinder
void Simulation::inject_particles(const int timestep) {
    static ThreadRandomGenerator randGenSpace;
    static ThreadRandomGenerator randGenPulse;
    const double3 cellSize = domain.cell_size();
    const double cellVolume = cellSize.x() * cellSize.y() * cellSize.z();
    const double dt = parameters.get_double("Dt");
    const int NumPartPerCell = parameters.get_int("NumPartPerCell");
    double3 center;
    center.x() = 0.5 * cellSize.x() * domain.num_cells(X);
    center.y() = 0.5 * cellSize.y() * domain.num_cells(Y);
    center.z() = 0.5 * cellSize.z() * domain.num_cells(Z);
    double rz = 0.5 * cellSize.z() * domain.num_cells(Z);
    double rr = 30 * cellSize.x();
    int ParticlesPerStep = dt * M_PI * rr * rr * (2 * rz) * NumPartPerCell /
                           (parameters.get_double("Tau") * cellVolume);
    for (auto &sp : species) {
        randGenSpace.SetRandSeed(100 + timestep);
        sp.injectionEnergy = sp.add_uniform_cilinderZ(
            ParticlesPerStep,
            double3(sp.temperature, sp.temperature, sp.temperature), center, rr,
            rz, randGenSpace, randGenPulse);
    }
}


// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each field stored in 1D array with 4d index x,y,z,d)
void Simulation::make_all() {
    bounds.setBounds(Bounds::BoundValues(BoundType::PERIODIC, BoundType::PERIODIC,
                               BoundType::PERIODIC),
                     Bounds::BoundValues(BoundType::PERIODIC, BoundType::PERIODIC,
                               BoundType::PERIODIC));
    domain.setDomain(parameters, bounds);
    const double dt  = parameters.get_double("Dt");
    const int startTimeStep = parameters.get_int("StartTimeStep");
    const int lastTimestep = parameters.get_int("LastTimestep");
    
    // allocate memory for fields
    init(mesh);
    
    
    Writer writer(mesh, species, domain, parameters);

    writer.output(0.0, parameters, startTimeStep);

    Timer globalTimer("globalFunctions.time");
    double energyP, energyPn;
    for (auto timestep = startTimeStep + 1; timestep <= lastTimestep;
         ++timestep) {
        globalTimer.start("Total");

        mesh.prepare();
        
        inject_particles(timestep);

        for (auto &sp : species) {
            sp.delete_bounds();
            sp.prepare(timestep);   // save start coord for esirkepov current
        }

        globalTimer.start("densityCalc");
        for (auto &sp : species) {
            sp.density_on_grid_update(); // calculate dendity field
        }
        globalTimer.finish("densityCalc");

        collect_charge_density(mesh.chargeDensityOld, species);

        globalTimer.start("particles1");
        for (auto &sp : species) {
            //sp.move(1. * Dt);   // +++ x_{n+1/2} -> x_{n+1}
            
            sp.move_and_calc_current(0.5 * dt);   //  +++ x_n -> x_{n+1/2}

            std::cout << "Moved " << sp.get_total_num_of_particles() << " "
                      << sp.name() << "\n";
            
            sp.update_cells(domain);
            // +++ get J(x_{n+1/2},v_n)_predict
            
            sp.predict_current(mesh.fieldB, mesh.fieldJp, domain, dt);
            
            // Stencil LmatX in special format. Very slow function
            // +++ get Lgg'(x_{n+1/2})
            sp.get_L(mesh, domain, dt);
        }
        globalTimer.finish("particles1");
        
        // mergy bounds for periodic 
        mesh.glue_Lmat_bound();

        globalTimer.start("stencilLmat");
        // convert LmatX to Eigen sparce matrix format Lmat
        mesh.stencil_Lmat(domain);
        globalTimer.finish("stencilLmat");

        mesh.make_periodic_border_with_add(mesh.fieldJp);

        globalTimer.start("FieldsPredict");
        // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})). mesh consist En,
        // En+1_predict
        mesh.predictE(dt); // solve system of linear equations using Lmat for find fieldE

        globalTimer.finish("FieldsPredict");

        globalTimer.start("particles2");
        energyP = energyPn = 0.;
        for (auto &sp : species) {
            energyP += sp.get_kinetic_energy();
            // +++ get v'_{n+1} from v_{n} and E'_{n+1}
            sp.predict_velocity(mesh, domain, dt); // calc new particles velocity using new fieldE
            // +++ x_{n+1/2} -> x_{n+1}
            //sp.move(0.5 * Dt);
            sp.move_and_calc_current(0.5 * dt);
            sp.update_cells(domain);
            mesh.make_periodic_border_with_add(sp.currentOnGrid);
            sp.currentOnGrid.data() *= 0.5;
        }
        globalTimer.finish("particles2");

        collect_current(mesh.fieldJe, species);

        globalTimer.start("FieldsCorr");
        // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final
        mesh.correctE(dt); // solve simple systeomof linear equations for correct fieldE
        globalTimer.finish("FieldsCorr");

        globalTimer.start("particles3");
        for (auto &sp : species) {
            sp.correctv(mesh, domain, dt); // correct particles velocity for save energy
            energyPn += sp.get_kinetic_energy();   // / sp.mpw(0);
        }
        for (auto &sp : species) {
            sp.density_on_grid_update();
        }
        globalTimer.finish("particles3");

        globalTimer.start("computeB");
        // calculate fieldB
        mesh.computeB(mesh.fieldE, mesh.fieldEn, mesh.fieldB, parameters.get_double("Dt"));
        globalTimer.finish("computeB");


        // later output data and check conservation layws
        collect_charge_density(mesh.chargeDensity, species);
        std::cout << mesh.chargeDensity.data().norm()
                  << " norm mesh.chargeDensity \n";

        auto divJ = mesh.divE * mesh.fieldJe.data();

        auto delta =
            (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) +
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
