// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <fstream>

#include "simulation.h"

// Main function simply hands off control to the Simulation class
int main(int argc, char **argv) {
    ParametersMap parameters(load_parameters("./SysParams.cfg"));
    ParametersMap outputParameters(load_parameters("./Diagnostics.cfg"));

    std::ifstream file("particles_config.json");
    if (!file.is_open()) {
        std::cerr << "Cannot open JSON file: " << "particles_config.json"
                  << std::endl;
        return 0;
    }

    nlohmann::json particles_config;
    file >> particles_config;

    auto simulation = build_simulation(parameters, particles_config,
                                       outputParameters, argc, argv);

    simulation->init();
    simulation->calculate();
    simulation->finalize();

    return 0;
}
