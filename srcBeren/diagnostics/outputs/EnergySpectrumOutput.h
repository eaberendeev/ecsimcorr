#pragma once

#include <string>
#include <vector>

#include "IDiagnosticOutput.h"
#include "ParticlesArray.h"

// === Energy spectrum text file per species ===
class EnergySpectrumOutput : public IDiagnosticOutput {
    const Species& species_;
    std::vector<std::string> species_names_;
    int interval_;

   public:
    EnergySpectrumOutput(const Species& species, const std::vector<std::string>& species_names, int interval);
    void output(int timestep, const Diagnostics&) override;
};
