// Author: Evgeny Berendeev

#include <map>
#include <string>

#include "Damping.h"
#include "Diagnostic.h"
#include "Mesh.h"
#include "Particles.h"
#include "Read.h"
#include "Vec.h"
#include "World.h"


void get_charge_density(Mesh &mesh, Field3d &field,
                        const std::vector<ParticlesArray> &species);
void get_current(Field3d &field, const std::vector<ParticlesArray> &species);

int main(int argc, char **argv) {
    Region regionGlob;

    Region region = regionGlob;
    /// Make modeling area
    World world(regionGlob, region);

    ///	Make mesh
    Mesh mesh(world);

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
            sp.prepare(timestep);   // save start coord for esirkepov current
        }

        globalTimer.start("densityCalc");
        for (auto &sp : species) {
            sp.density_on_grid_update();
        }
        globalTimer.finish("densityCalc");

        get_charge_density(mesh, mesh.chargeDensityOld, species);

        globalTimer.start("particles1");
        for (auto &sp : species) {
            // sp.move(0.5*Dt); //  +++ x_n -> x_{n+1/2}
            sp.move_and_calc_current(0.5 * Dt);   //  +++ x_n -> x_{n+1/2}
            sp.predict_current(
                mesh.fieldB,
                mesh.fieldJp);   // +++ get J(x_{n+1/2},v_n)_predict
            sp.update_cells();

            sp.get_L(mesh);   // +++ get Lgg'(x_{n+1/2})
        }
        globalTimer.finish("particles1");

        mesh.glue_Lmat_bound();

        globalTimer.start("stencilLmat");
        mesh.stencil_Lmat();
        globalTimer.finish("stencilLmat");

        mesh.make_periodic_border_with_add(mesh.fieldJp);

        globalTimer.start("FieldsPredict");
        mesh.predictE();   // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})).
                           // mesh consist En, En+1_predict
        globalTimer.finish("FieldsPredict");

        globalTimer.start("particles2");
        energyP = energyPn = 0.;
        for (auto &sp : species) {
            energyP += sp.get_kinetic_energy();
            sp.predict_velocity(
                mesh);   // +++ get v'_{n+1} from v_{n} and E'_{n+1}

            // sp.move(0.5*Dt); // +++ x_{n+1/2} -> x_{n+1}
            sp.move_and_calc_current(0.5 * Dt);   //  +++ x_n -> x_{n+1/2}

            sp.update_cells();
            mesh.make_periodic_border_with_add(sp.currentOnGrid);
            sp.currentOnGrid.data() *= 0.5;
        }
        globalTimer.finish("particles2");

        get_current(mesh.fieldJe, species);

        globalTimer.start("FieldsCorr");
        mesh.correctE();   // ---- get E_{n+1} from E_n and J_e. mesh En changed
                           // to En+1_final
        globalTimer.finish("FieldsCorr");

        globalTimer.start("particles3");
        for (auto &sp : species) {
            sp.correctv(mesh);
            energyPn += sp.get_kinetic_energy();   // / sp.mpw(0);
        }
        for (auto &sp : species) {
            sp.update_cells();
            sp.density_on_grid_update();
        }
        globalTimer.finish("particles3");

        globalTimer.start("computeB");
        mesh.computeB();
        globalTimer.finish("computeB");

        get_charge_density(mesh, mesh.chargeDensity, species);
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

        damping_fields(mesh.fieldEn, mesh.fieldB, region);

        mesh.fieldB.data() += mesh.fieldBInit.data();

        /// END DAMPING ///////
        globalTimer.finish("Total");
        globalTimer.write(timestep);
        for (auto &sp : species) {
            sp.write_particles_to_recovery(timestep);
        }
        mesh.write_fields_to_recovery(timestep);
    }

    return 0;
}

void get_charge_density(Mesh &mesh, Field3d &field,
                        const std::vector<ParticlesArray> &species) {
    field.clear();
    for (const auto &sp : species) {
#pragma omp parallel for collapse(3)
        for (auto i = 0; i < sp.densityOnGrid.size().x(); i++) {
            for (auto j = 0; j < sp.densityOnGrid.size().y(); j++) {
                for (auto k = 0; k < sp.densityOnGrid.size().z(); k++) {
                    field(mesh.sind(i, j, k)) += sp.densityOnGrid(i, j, k);
                }
            }
        }
    }
}

void get_current(Field3d &field, const std::vector<ParticlesArray> &species) {
    field.clear();
    for (const auto &sp : species) {
        field.data() += sp.currentOnGrid.data();
    }
}
