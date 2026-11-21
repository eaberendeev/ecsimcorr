#pragma once

#include <string>
#include <vector>

#include "diagnostics/DiagnosticMeasurements.h"
#include "diagnostics/IDiagnosticOutput.h"

// Velocity distribution f(a, b) of species in a spatial box, one text file per
// species and output step: Particles/<name>/Distribution/<basename>_<step>.txt.
class DistributionOutput : public IDiagnosticOutput {
   public:
    DistributionOutput(const std::vector<std::string>& species_names, const DistributionRequest& request, int interval,
                       std::string basename);
    void request_measurements(int timestep, DiagnosticMeasurements& measured) const override;
    void output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) override;

   private:
    std::vector<std::string> species_names_;
    DistributionRequest request_;
    int interval_;
    std::string basename_;
};
