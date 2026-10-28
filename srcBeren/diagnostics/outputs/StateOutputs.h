#pragma once

#include <string>
#include <vector>

#include "DiagnosticOutputSettings.h"
#include "IDiagnosticOutput.h"
#include "ParticlesArray.h"
#include "containers.h"

// === Checkpoint / recovery ===
class RecoveryOutput : public IDiagnosticOutput {
    const Field3d& fieldEn_;
    const Field3d& fieldBn_;
    const Species& species_;
    int interval_;

   public:
    RecoveryOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Species& species, int interval);
    void output(int timestep, const Diagnostics&) override;
};

// === 3D field dump at specific timesteps ===
class Fields3DOutput : public IDiagnosticOutput {
    std::vector<OutputField> fields_;
    std::vector<int> timesteps_;

   public:
    Fields3DOutput(const std::vector<OutputField>& fields, const std::vector<int>& timesteps);
    void output(int timestep, const Diagnostics&) override;
};

// === 3D all-field + per-species dump at specific timesteps ===
// Writes fields (E/B/En/Bn/E_ex/B_init) + per-species density/current as 3D binary files
class Full3DAllOutput : public IDiagnosticOutput {
    const Field3d& fieldEn_;
    const Field3d& fieldBn_;
    const Field3d& fieldE_external_;
    const Field3d& fieldBInit_;
    const Species& species_;
    std::vector<int> timesteps_;

   public:
    Full3DAllOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                    const Field3d& fieldBInit, const Species& species, const std::vector<int>& timesteps);
    void output(int timestep, const Diagnostics&) override;
};
