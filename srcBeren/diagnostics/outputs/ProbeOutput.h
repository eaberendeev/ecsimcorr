#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "IDiagnosticOutput.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"

// === Point probes: interpolates all fields at given (x,y,z) points, writes probes.txt ===
class ProbeOutput : public IDiagnosticOutput {
    const Domain& domain_;
    const Field3d& fieldEn_;
    const Field3d& fieldBn_;
    const Field3d& fieldE_external_;
    const Field3d& fieldBInit_;
    const Species& species_;
    std::vector<std::string> species_names_;
    std::vector<Vector3R> points_;
    int interval_;
    std::ofstream fProbes_;
    bool header_written_;

   public:
    ProbeOutput(const Domain& domain, const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                const Field3d& fieldBInit, const Species& species, const std::vector<std::string>& species_names,
                const std::vector<Vector3R>& points, int interval);
    void output(int timestep, const Diagnostics&) override;
    void finalize() override;
};
