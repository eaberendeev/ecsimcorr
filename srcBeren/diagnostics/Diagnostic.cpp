#include "Diagnostic.h"

Diagnostics::Diagnostics(const Species& species) {
    for (const auto& kv : species) particleNames.push_back(kv.second->name());
}

double Diagnostics::value_or_zero(const std::string& key) const {
    const auto it = values.find(key);
    return it == values.end() ? 0.0 : it->second;
}

void Diagnostics::addValue(const std::string& key, double value) {
    if (values.find(key) == values.end()) {
        valueOrder.push_back(key);
    }
    values[key] = value;
}

void Diagnostics::addBoundary(const std::string& key, double value) {
    if (values.find(key) == values.end()) {
        boundaryOrder.push_back(key);
    }
    values[key] = value;
}
