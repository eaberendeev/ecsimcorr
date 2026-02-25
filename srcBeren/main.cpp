// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <fstream>

#include "simulation.h"

// Main function simply hands off control to the Simulation class
int main(int argc, char** argv) {
    std::ifstream system_config_file("system_config.json");
    std::ifstream particles_config_file("particles_config.json");

    if (!system_config_file.is_open() || !particles_config_file.is_open()) {
        std::cerr << "Cannot open JSON file: "
                  << "system_config.json / particles_config.json" << std::endl;
        return 0;
    }

    nlohmann::json system_config;
    nlohmann::json particles_config;
    system_config_file >> system_config;
    particles_config_file >> particles_config;

    auto simulation = build_simulation(system_config,
                                       particles_config, argc, argv);

    simulation->init();
    simulation->calculate();
    simulation->finalize();

    return 0;
}
