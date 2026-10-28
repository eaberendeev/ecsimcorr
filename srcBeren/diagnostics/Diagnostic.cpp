#include "Diagnostic.h"

Diagnostics::Diagnostics(const Species& species) {
    for (const auto& kv : species) particleNames.push_back(kv.second->name());
}

double Diagnostics::value_or_zero(const std::string& key) const {
    const auto it = energy.find(key);
    return it == energy.end() ? 0.0 : it->second;
}

void Diagnostics::addEnergy(const std::string& key, double value) {
    if (energy.find(key) == energy.end()) {
        energyOrder.push_back(key);
    }
    energy[key] = value;
}

void Diagnostics::addBoundary(const std::string& key, double value) {
    if (energy.find(key) == energy.end()) {
        boundaryOrder.push_back(key);
    }
    energy[key] = value;
}
