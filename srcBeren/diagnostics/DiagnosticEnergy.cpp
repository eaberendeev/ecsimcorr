#include <assert.h>

#include "Diagnostic.h"

void Diagnostics::addEnergy(const std::string &key, double value) {
    if (energy.find(key) == energy.end()) {
        energyOrder.push_back(key);
    }
    energy[key] = value;
}

void Diagnostics::addBoundary(const std::string &key, double value) {
    if (energy.find(key) == energy.end()) {
        boundaryOrder.push_back(key);
    }
    energy[key] = value;
}
