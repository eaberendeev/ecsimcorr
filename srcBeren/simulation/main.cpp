// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <fstream>

#include "simulation_implicit.h"
//#include "simulation_ecsim_corr.h"

// Main function simply hands off control to the Simulation class
int main(int argc, char **argv) {
    ParametersMap parameters(load_parameters("./SysParams.cfg"));
    std::vector<ParametersMap> speciesParameters =
    load_vector_parameters("./PartParams.cfg", "Particles");
    ParametersMap outputParameters(load_parameters("./Diagnostics.cfg"));

    SimulationImplicit simulation(parameters, speciesParameters,
                                  outputParameters, argc, argv);
    //SimulationEcsimCorr simulation(argc, argv);
    simulation.init();
    simulation.make_all();

    return 0;
}
