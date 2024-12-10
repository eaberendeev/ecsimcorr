// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <fstream>

#include "test_track_in_fieldB.h"

void set_test_parameters(ParametersMap& parameters,
                         std::vector<ParametersMap>& speciesParameters,
                         ParametersMap& outputParameters) {
    //  parameters
    parameters.set("SystemParams", {"SystemParams"});
    parameters.set("Dx", {"0.5"});
    parameters.set("Dy", {"0.5"});
    parameters.set("Dz", {"0.5"});
    parameters.set("Dt", {"5."});
    parameters.set("NumCellsX_glob", {"5"});
    parameters.set("NumCellsY_glob", {"5"});
    parameters.set("NumCellsZ_glob", {"5"});
    parameters.set("DampCellsX_glob", {"0", "0"});
    parameters.set("DampCellsY_glob", {"0", "0"});
    parameters.set("DampCellsZ_glob", {"0", "0"});

    parameters.set("BoundTypeX", {"PERIODIC", "PERIODIC"});
    parameters.set("BoundTypeY", {"PERIODIC", "PERIODIC"});
    parameters.set("BoundTypeZ", {"PERIODIC", "PERIODIC"});
    parameters.set("BCoil", {"0"});
    parameters.set("MC2", {"512"});
    parameters.set("StartFromTime", {"0"});
    parameters.set("Collider", {"None"});
    parameters.set("DampingType", {"None"});
    parameters.set("NumPartPerCell", {"50"});

    parameters.set("BUniform", {"0", "0", "0"});
    parameters.set("RecoveryInterval", {"800"});
    parameters.set("TimeStepDelayDiag1D", {"20001"});
    parameters.set("TimeStepDelayDiag2D", {"8000"});
    parameters.set("n0", {"10000000000000.0"});
    parameters.set("k_particles_reservation", {"-1.0"});
    parameters.set("StartTimeStep", {"0"});

    parameters.set("LastTimestep", {"100"});

    // speciesParameters

    ParametersMap electronParameters;
    electronParameters.set("Particles", {"Electrons"});
    electronParameters.set("Charge", {"-1.0"});
    electronParameters.set("Density", {"1.0"});
    electronParameters.set("Velocity", {"0.0"});
    electronParameters.set("Mass", {"1.0"});
    electronParameters.set("Temperature", {"0", "0.0", "0.0"});
    electronParameters.set("DistType", {"None"});
    electronParameters.set("DistSpace", {"None"});
    electronParameters.set("DistPulse", {"Gauss"});
    speciesParameters.emplace_back(electronParameters);

    // outputParameters
    outputParameters.set("Diagnostics", {"Diagnostics"});
}


void particle_trajectory_test_in_constant_fields(ParametersMap& parameters,
            std::vector<ParametersMap>& speciesParameters,
            ParametersMap& outputParameters, int argc, char** argv) {
    outputParameters.set("TestType", {"Bz"});
    outputParameters.set("MoverType", {"ecsim"});
    double MaxTime = 200;
    Bounds bounds;
    Domain domain;
    bounds.setBounds(parameters);
    domain.setDomain(parameters, bounds);
    Field3d B(domain.size(), 3);
    Field3d E(domain.size(), 3);
    /// set fields set particles
    for (int i = 0; i < domain.size().x(); i++) {
        for (int j = 0; j < domain.size().y(); j++) {
            for (int k = 0; k < domain.size().z(); k++) {
                B(i, j, k, Axis::Z) =
                    0.2 *
                    (1 - 0.8 * (domain.cell_size(Axis::Y) * (j - 0.5) - 1.25));
            }
        }
    }
    E.set_zero();

    Particle ptest(1.25, 1.25 - 0.1, 1.25, 0.02, 0.0, 0.0);
    ptest.id = 0;
    std::vector<std::pair<std::string, Particle>> pairs;
    pairs.push_back(std::make_pair("Electrons",ptest));
    std::vector<std::string> dt{"0.5", "1.5", "3"};

    for(auto& d:dt){

        parameters.set("Dt", {d});
        double ddt = stod(d);
        parameters.set("LastTimestep", {to_string(MaxTime/ddt,4)});

        ConstantFieldParticleTrajectorySimulator simulation(
            parameters, speciesParameters, outputParameters, argc, argv);
        simulation.init();
        simulation.init_fields(E, B);
        simulation.init_particles(pairs);
        simulation.init_diagnostic();
        simulation.make_all();
    }
}

int main(int argc, char** argv) {
    ParametersMap parameters;
    std::vector<ParametersMap> speciesParameters;
    ParametersMap outputParameters;
    set_test_parameters(parameters, speciesParameters, outputParameters);

    particle_trajectory_test_in_constant_fields(parameters, speciesParameters,
                                                outputParameters, argc, argv);

    return 0;
}
