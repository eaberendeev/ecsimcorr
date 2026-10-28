#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "IDiagnosticOutput.h"
#include "containers.h"

// === Point probes: writes prepared fields and densities at given (x,y,z) points ===
class ProbeOutput : public IDiagnosticOutput {
    std::vector<std::string> species_names_;
    std::vector<Vector3R> points_;
    int interval_;
    std::ofstream fProbes_;
    bool header_written_;

   public:
    ProbeOutput(const std::vector<std::string>& species_names, const std::vector<Vector3R>& points, int interval);
    void request_measurements(int timestep, DiagnosticMeasurements& measured) const override;
    void output(int timestep, const Diagnostics&, const DiagnosticMeasurements&) override;
    void finalize() override;
};
