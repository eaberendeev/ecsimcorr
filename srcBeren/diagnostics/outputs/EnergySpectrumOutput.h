#pragma once

#include <string>
#include <vector>

#include "IDiagnosticOutput.h"

// === Energy spectrum text file per species ===
class EnergySpectrumOutput : public IDiagnosticOutput {
    std::vector<std::string> species_names_;
    int interval_;

   public:
    EnergySpectrumOutput(const std::vector<std::string>& species_names, int interval);
    void request_measurements(int timestep, DiagnosticMeasurements& measured) const override;
    void output(int timestep, const Diagnostics&, const DiagnosticMeasurements&) override;
};
