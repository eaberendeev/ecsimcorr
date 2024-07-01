// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <iostream>

#include "recovery.h"
#include "simulation.h"
#include "util.h"
#include "service.h"

void Simulation::output_all(const Mesh& mesh,
                            const std::vector<ParticlesArray>& species,
                            int timestep) {
    output_energy(timestep);
        // Recovery
        for (auto& sp : species) {
        write_particles_to_recovery(sp, timestep, parameters.get_int("RecoveryInterval"));
    }
}

void Simulation::output_energy(const int timestep){
    std::stringstream ss;
    static FILE* fDiagEnergies = fopen("Energies.dat", "w");
    static std::map<std::string,double> energy;

    for (auto& sp : species) {
        energy[sp.name() + "Area"] = sp.get_kinetic_energy();
        energy[sp.name() + "Z"] = sp.get_kinetic_energy(Z);
        energy[sp.name() + "XY"] = sp.get_kinetic_energy(X, Y);
  }

    if (timestep == 0) {
        ss << "Time ";
        for (auto it = energy.begin(); it != energy.end();
             ++it) {
            ss << "Energy_" << it->first << " ";
        }
        ss << "\n";
    }

    ss << timestep * parameters.get_double("Dt") << " ";

    for (auto it = energy.begin(); it != energy.end(); ++it) {
        double energyP = it->second;
        ss << energyP << " ";
    }
    ss << "\n";
    fprintf(fDiagEnergies, "%s", (ss.str()).c_str());
    std::cout << ss.str();
    if( timestep % parameters.get_int("TimeStepDelayDiag1D") == 0){
      fflush(fDiagEnergies);
    }
}
void Simulation::make_folders() const {
    create_directory(".//Anime");
    create_directory(".//Recovery");
    create_directory(".//Recovery//Particles");
    create_directory(".//Performance");
    create_directory(".//Particles");
    for (const auto& sp : species) {
        create_directory(".//Particles//" + sp.name());
        create_directory(".//Recovery//Particles//" + sp.name());
    }
    std::cerr << "Folders for output has been created: SUCCESS\n";
}
