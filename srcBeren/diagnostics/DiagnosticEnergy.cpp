#include <assert.h>

#include "Diagnostic.h"

void Diagnostics::addEnergy(const std::string &key, double value) {
    if (energy.find(key) == energy.end()) {
        energyOrder.push_back(key);
    }
    energy[key] = value;
}

void Diagnostics::write_energy(const nlohmann::json &system_config, int timestep) {
    std::stringstream ss;

    static bool writeHeader = false;
    if (!writeHeader) {
        ss << "Time ";
        for (const auto &key : energyOrder) {
            ss << key << " ";
        }
        writeHeader = true;
        ss << "\n";
    }
    ss << timestep*get_checked<double>(system_config, "Dt") << " ";
    for (const auto &key : energyOrder) {
        ss << energy[key] << " ";
    }
    ss << "\n";
    std::cout << ss.str();
    fEnergy << ss.str();
    fEnergy.flush();
}