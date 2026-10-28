#pragma once

class Diagnostics;

// === Base class for all configurable diagnostic outputs ===
// Pattern: virtual method + factory, like IDistribution in particles_distribution_collection.h
class IDiagnosticOutput {
   public:
    virtual ~IDiagnosticOutput() = default;
    virtual void output(int timestep, const Diagnostics& energy_diag) = 0;
    // Called once after the simulation loop: flushes residual averaged windows
    // and closes any files. Must be idempotent.
    virtual void finalize() {
    }
};
