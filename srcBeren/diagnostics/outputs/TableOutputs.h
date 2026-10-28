#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "Diagnostic.h"
#include "IDiagnosticOutput.h"

// === Shared table output (energy.txt / boundary.txt) ===
class TableOutput : public IDiagnosticOutput {
   public:
    void output(int timestep, const Diagnostics& energy_diag, const DiagnosticMeasurements&) override;
    void finalize() override;

   protected:
    TableOutput(const std::string& filename, double dt, int interval);
    virtual const std::vector<std::string>& column_order(const Diagnostics& diag) const = 0;

   private:
    std::ofstream file_;
    std::string filename_;
    double dt_;
    int interval_;
    bool header_written_;
};

// === Energy balance (energy.txt) ===
class EnergyBalanceOutput : public TableOutput {
   public:
    EnergyBalanceOutput(double dt, int interval) : TableOutput("energy.txt", dt, interval) {
    }

   protected:
    const std::vector<std::string>& column_order(const Diagnostics& diag) const override {
        return diag.energyOrder;
    }
};

// === Boundary stats (boundary.txt) ===
class BoundaryOutput : public TableOutput {
   public:
    BoundaryOutput(double dt, int interval) : TableOutput("boundary.txt", dt, interval) {
    }

   protected:
    const std::vector<std::string>& column_order(const Diagnostics& diag) const override {
        return diag.boundaryOrder;
    }
};

// === Console summary (per-step print) ===
class ConsoleOutput : public IDiagnosticOutput {
    double dt_;

   public:
    ConsoleOutput(double dt);
    void output(int timestep, const Diagnostics& energy_diag, const DiagnosticMeasurements&) override;
};

// Full-grid norms of (rho^{n+1}-rho^n)/dt + div J for selected Corr steps.
class ChargeContinuityOutput : public IDiagnosticOutput {
    std::ofstream file_;
    double dt_;
    int interval_;

   public:
    ChargeContinuityOutput(double dt, int interval);
    bool checks_charge_continuity(int step) const override {
        return step > 0 && step % interval_ == 0;
    }
    void output(int timestep, const Diagnostics& values, const DiagnosticMeasurements&) override;
    void finalize() override {
        if (file_.is_open())
            file_.close();
    }
};
