#ifndef DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "ParticlesArray.h"

struct ChargeContinuityResidual {
    double l2;
    double linf;
};

class Diagnostics {
   public:
    explicit Diagnostics(const Species& species);

    // Output preserves the legacy zero for a missing value without inserting
    // a new entry into the measured diagnostic state.
    double value_or_zero(const std::string& key) const;

    // Energy tracking API
    void addEnergy(const std::string& key, double value);
    void addBoundary(const std::string& key, double value);

    // Data for outputs to read
    std::map<std::string, double> energy;
    std::vector<std::string> energyOrder;
    std::vector<std::string> boundaryOrder;
    std::vector<std::string> particleNames;
    // Present only on selected completed ECSIM-Corr steps; never persisted.
    std::optional<ChargeContinuityResidual> charge_continuity;
};

#endif
