#pragma once

class Diagnostics;
class DiagnosticMeasurements;

// === Base class for all configurable diagnostic outputs ===
// Pattern: virtual method + factory, like IDistribution in particles_distribution_collection.h
class IDiagnosticOutput {
   public:
    virtual ~IDiagnosticOutput() = default;
    // Register only measurements needed by this output at the current step.
    virtual void request_measurements(int, DiagnosticMeasurements&) const {
    }
    // Continuity needs both endpoints of the selected step, before output time.
    virtual bool checks_charge_continuity(int) const {
        return false;
    }
    virtual void output(int timestep, const Diagnostics& values, const DiagnosticMeasurements& measured) = 0;
    // Called once after the simulation loop: flushes residual averaged windows
    // and closes any files. Must be idempotent.
    virtual void finalize() {
    }
};
