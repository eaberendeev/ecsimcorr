// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <fstream>

#include "simulation.h"

// Main function simply hands off control to the Simulation class
int main(int argc, char **argv) {
    ParametersMap parameters(load_parameters("./SysParams.cfg"));
    std::vector<ParametersMap> speciesParameters =
    load_vector_parameters("./PartParams.cfg", "Particles");
    ParametersMap outputParameters(load_parameters("./Diagnostics.cfg"));

    auto simulation =
        build_simulation(parameters, speciesParameters, outputParameters, argc, argv);

    simulation->init();
    simulation->calculate();
    simulation->finalize();

    return 0;
}
