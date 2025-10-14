// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_symplectic.h"

#include <iostream>
#include <map>
#include <string>

#include "Coil.h"
#include "Damping.h"
#include "Diagnostic.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "ParticlesDiagnostic.h"
#include "Read.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "recovery.h"

// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each field stored in 1D array with 4d index x,y,z,d)
void SimulationSymplectic::make_step([[maybe_unused]] const int timestep) {
    const double dt = parameters.get_double("Dt");
    globalTimer.start("Total");
    double ek = 0, ek1 = 0, ek2 = 0;
    globalTimer.start("densityCalc");
    for (auto &sp : species) {
        if (sp->name() == "Electrons" && timestep == 1) {
            Particle ptest(2.195, 2.19 , 2.1, 0.1, 0.2, 0.0);
            sp->add_particle(ptest);
        }
        sp->density_on_grid_update(SHAPE);   // calculate dendity field
        ek += sp->get_kinetic_energy();
    }
    globalTimer.finish("densityCalc");

    collect_charge_density(mesh.chargeDensityOld);
    mesh.apply_density_boundaries(mesh.chargeDensityOld, domain);
    fieldBFull = fieldB + fieldBInit;

    move_particles(fieldBFull, fieldJ);

    mesh.apply_density_boundaries(fieldJ, domain);
    std::cout << "fieldJ.norm() = " << fieldJ.norm() << "\n";
    std::cout << "fieldE.norm() = " << fieldE.norm() << "\n";

    Field3d fieldE0 = fieldE;
    double energy0 = mesh.calc_energy_field(fieldE0);

    update_E();
    for (auto &sp : species) {
        ek1 += sp->get_kinetic_energy();

        update_velocity(sp, dt);
        ek2 += sp->get_kinetic_energy();
    }
    for (auto &sp : species) {
        sp->update_cells(domain);
    }
    update_B();

    Field3d ff = 0.5 * (fieldE0 + fieldE);
    const double energyJeEn = calc_JE(ff, fieldJ, bounds);
    double energy = mesh.calc_energy_field(fieldE);

    std::cout << "Energy JE: " << dt*energyJeEn << " "
              << mesh.calc_energy_field(fieldB) + energy - energy0 << " " << ek2 - ek << std::endl;
    std::cout << "energyK = " << ek << " " << ek1 <<  " " << ek2  << "\n";
    std::cout << "v2  = " << sqrt(ek2*2.0) << "\n";
    for (auto &sp : species) {
        sp->density_on_grid_update(SHAPE);
    }
    std::cout << "density_on_grid_update " << "\n";
    // later output data and check conservation layws
    collect_charge_density(mesh.chargeDensity);
    mesh.apply_density_boundaries(mesh.chargeDensity, domain);

    auto divJ = mesh.divE * fieldJ.data();

    auto delta =
        (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) +
        divJ;
    std::cout << delta.norm() << " norm drho / Dt - divJ \n";
    globalTimer.finish("Total");
    for (int i = 0; i < mesh.chargeDensity.size().x()) {
        for(int j = 0; j < mesh.chargeDensity.size().y()){
            for (int k = 0; k < mesh.chargeDensity.size().z()) {
                if(mesh.chargeDensity(i,j,k,0) != 0)
                    std::cout << i << " " << j << " " << k << " "
                              << mesh.chargeDensity(i, j, k, 0) << "\n";
            }
         }
    }
    for (int i = 0; i < mesh.chargeDensity.size().x()) {
        for (int j = 0; j < mesh.chargeDensity.size().y()) {
            for (int k = 0; k < mesh.chargeDensity.size().z()) {
                if (mesh.chargeDensity(i, j, k, 0) != 0)
                    std::cout << i << " " << j << " " << k << " "
                              << mesh.chargeDensityOld(i, j, k, 0) << "\n";
            }
        }
    }
    for (int i = 0; i < fieldJ.size().x()) {
        for (int j = 0; j < fieldJ.size().y()) {
            for (int k = 0; k < fieldJ.size().z()) {
                if (fieldJ(i, j, k, X) != 0)
                    std::cout << i << " " << j << " " << k << " "
                              << fieldJ(i, j, k, X) << "\n";
            }
        }
    }
    for (int i = 0; i < fieldJ.size().x()) {
        for (int j = 0; j < fieldJ.size().y()) {
            for (int k = 0; k < fieldJ.size().z()) {
                if (fieldJ(i, j, k, Y) != 0)
                    std::cout << i << " " << j << " " << k << " "
                              << fieldJ(i, j, k, Y) << "\n";
            }
        }
    }
     if (timestep > 1)
         exit(0);
}

void SimulationSymplectic::init_fields() {
    fieldJ.resize(domain.size(), 3);
    fieldE.resize(domain.size(), 3);
    fieldB.resize(domain.size(), 3);
    fieldBInit.resize(domain.size(), 3);
    fieldBFull.resize(domain.size(), 3);

    fieldJ.setZero();

    fieldE.setZero();
    fieldB.setZero();

    if (parameters.get_int("StartFromTime") > 0) {
        read_fields_from_recovery(fieldE, fieldB);
    } else {
        mesh.set_uniform_field(fieldE, 0.0, 0.0, 0.0);
    }
    mesh.set_uniform_field(fieldBInit, parameters.get_double("BUniform", 0),
                           parameters.get_double("BUniform", 1),
                           parameters.get_double("BUniform", 2));
    set_coils(fieldBInit, domain, parameters);

}

void SimulationSymplectic::prepare_step(const int timestep) {
    inject_particles(timestep, domain);
    //damping_fields(fieldE, fieldB, domain,
    //               parameters);
    fieldJ.setZero();

    for (auto &sp : species) {
        sp->prepare();   // save start coord for esirkepov current
        std::cout << sp->get_total_num_of_particles() << " size \n";
    }
}

void SimulationSymplectic::make_diagnostic(const int timestep) {
    static Diagnostics diagnostic(outputParameters, domain, species);

    for (auto &sp : species) {
        write_particles_to_recovery(sp, timestep,
                                    parameters.get_int("RecoveryInterval"));
    }
    write_fields_to_recovery(fieldE, fieldB, timestep,
                             parameters.get_int("RecoveryInterval"));
    // writer.output(0, parameters, timestep);
    diagnostic_energy(diagnostic);
    diagnostic.write_energy(parameters, timestep);
    const std::string pathToField = ".//Fields//Diag2D//";

    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    std::vector<std::pair<const Field3d &, std::string>> fields = {
        {fieldE, pathToField + "FieldE"},
        {fieldBFull, pathToField + "FieldB"}};
    diagnostic.output_fields2D(timestep, fields);
    for (auto &sp : species) {
        if (timestep % parameters.get_int("TimeStepDelayDiag2D") == 0) {
            const std::string spectrumPath =
                ".//Particles//" + sp->name() + "//";
            EnergySpectrum spectrum = sp->calculate_energy_spectrum();
            diagnostic.output_energy_spectrum(spectrum, timestep, spectrumPath);
        }

        const std::string pathToField =
            ".//Particles//" + sp->name() + "//Diag2D//";
        if(sp->is_neutral()){
            std::vector<std::pair<const Field3d &, std::string>> fields = {
                {sp->densityOnGrid, pathToField + "Density"}};
            diagnostic.output_fields2D(timestep, fields);
            continue;
        }
        Field3d pressureRR(domain.size(), 1);
        Field3d pressurePP(domain.size(), 1);
        Field3d pressureZZ(domain.size(), 1);
        Field3d pressureRP(domain.size(), 1);
        Field3d pressureRZ(domain.size(), 1);
        Field3d pressureZP(domain.size(), 1);
        // Использование:
        sp->calculate_pressure_component(pressureRR, RadialVelocity{}, RadialVelocity{});
        sp->calculate_pressure_component(pressurePP, PhiVelocity{}, PhiVelocity{});
        sp->calculate_pressure_component(pressureZZ, ZVelocity{}, ZVelocity{});
        sp->calculate_pressure_component(pressureRP, RadialVelocity{},
                                         PhiVelocity{});
        sp->calculate_pressure_component(pressureRZ, RadialVelocity{}, ZVelocity{});
        sp->calculate_pressure_component(pressureZP, ZVelocity{},
                                         PhiVelocity{});
        std::vector<std::pair<const Field3d &, std::string>> fields = {
            {sp->currentOnGrid, pathToField + "Current"},
            {sp->densityOnGrid, pathToField + "Density"},
            {pressureRR, pathToField + "Prr"},
            {pressurePP, pathToField + "Ppp"},
            {pressureZZ, pathToField + "Pzz"},
            {pressureRP, pathToField + "Prp"},
            {pressureRZ, pathToField + "Prz"},
            {pressureZP, pathToField + "Pzp"}};
        diagnostic.output_fields2D(timestep, fields);
    }
#ifdef SET_PARTICLE_IDS
    static ParticleTracker tracker(species, 1, "Tracking", "");
    tracker.track_particles(species, timestep);
#endif
}

void SimulationSymplectic::diagnostic_energy(Diagnostics &diagnostic) {
    // Статические переменные для сохранения состояния между вызовами
    static bool is_initialized = false;
    static double total_energy_initial = 0.0;
    static Field3d fieldE_initial;

    // 1. Вычисляем текущую кинетическую энергию частиц
    double kinetic_energy = 0;
    double kinetic_energy_inject = 0;
    double kinetic_energy_lost = 0;
    for (auto &sp : species) {
        // Сохраняем диагностические данные (оригинальные записи)
        diagnostic.addEnergy(sp->name(), sp->get_kinetic_energy());
        diagnostic.addEnergy(sp->name() + "Inject", sp->injectionEnergy);
        diagnostic.addEnergy(sp->name() + "LostZ", sp->lostEnergyZ);
        diagnostic.addEnergy(sp->name() + "LostXY", sp->lostEnergyXY);
        diagnostic.addEnergy(sp->name() + "Z", sp->get_kinetic_energy(Z));
        diagnostic.addEnergy(sp->name() + "XY", sp->get_kinetic_energy(X, Y));
        kinetic_energy += diagnostic.energy[sp->name()];
        kinetic_energy_inject += diagnostic.energy[sp->name() + "Inject"];
        kinetic_energy_lost += diagnostic.energy[sp->name() + "LostZ"] + diagnostic.energy[sp->name() + "LostXY"];
    }

    // 2. Вычисляем энергию полей
    fieldBFull.data() = fieldB.data() + fieldBInit.data();
    diagnostic.addEnergy("energyFieldE", mesh.calc_energy_field(fieldE));
    diagnostic.addEnergy("energyFieldB", mesh.calc_energy_field(fieldB));
    diagnostic.addEnergy("energyFieldBFull",
                         mesh.calc_energy_field(fieldBFull));
    double field_energy = diagnostic.energy["energyFieldE"] + diagnostic.energy["energyFieldB"];

    // 3. Рассчитываем полную энергию системы
    double total_energy = kinetic_energy + field_energy + kinetic_energy_lost; // - kinetic_energy_inject;

    // 4. Инициализация при первом вызове
    if (!is_initialized) {
        total_energy_initial = total_energy;
        is_initialized = true;
        fieldE_initial.resize(domain.size(), 3);
        fieldE_initial.setZero();
        // Записываем начальную энергию для последующего сравнения
        diagnostic.addEnergy("totalEnergyInitial", total_energy_initial);
    }

    // 5. Проверяем сохранение энергии
    double conservation_error = total_energy/total_energy_initial - 1;
    diagnostic.addEnergy("totalEnergy", total_energy);
    diagnostic.addEnergy("energyConservationError", conservation_error);

    // 6. Дополнительная диагностика (опционально)
    const double dt = parameters.get_double("Dt");
    diagnostic.addEnergy("JpEnergy", 0.5 * fieldJ.data().dot(
                                               fieldE.data() + fieldE_initial.data()));
    diagnostic.addEnergy("JpEnergyConservation",
                         diagnostic.energy["JpEnergy"] - field_energy);
    fieldE_initial = fieldE;
}