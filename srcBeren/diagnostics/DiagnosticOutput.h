#pragma once
#ifndef DIAGNOSTIC_OUTPUT_H_
#define DIAGNOSTIC_OUTPUT_H_

#include <array>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Diagnostic.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"

// === Base class for all configurable diagnostic outputs ===
// Pattern: virtual method + factory, like IDistribution in particles_distribution_collection.h
class IDiagnosticOutput {
   public:
    virtual ~IDiagnosticOutput() = default;
    virtual void output(int timestep, Diagnostics& energy_diag) = 0;
    // Called once after the simulation loop: flushes residual averaged windows
    // and closes any files. Must be idempotent.
    virtual void finalize() {
    }
};

// An averaging window stores only requested planes, including ghost nodes just
// like instantaneous files. All named fields share the output-call count.
class AveragedField2D {
   public:
    void set_planes(const std::vector<double>& px, const std::vector<double>& py, const std::vector<double>& pz,
                    const Domain& domain);
    void add(const std::string& key, const Field3d& field, const Field3d* extra = nullptr);
    void increment();
    int count() const;
    bool ready(int steps) const;
    void write(const std::string& key, const std::string& basename, int sequence, bool averaged = true) const;
    std::vector<std::string> keys() const;
    long long sequence() const;
    size_t stored_values() const;
    void reset();

   private:
    struct Plane {
        int axis;
        int position;
        Field3d sum;   // Size one in the normal direction; double precision until output.
    };
    std::array<std::vector<int>, 3> positions_;
    std::map<std::string, std::vector<Plane>> accum_;
    int count_ = 0;
    long long sequence_ = 0;
};

// === Shared table output (energy.txt / boundary.txt) ===
class TableOutput : public IDiagnosticOutput {
   public:
    void output(int timestep, Diagnostics& energy_diag) override;
    void finalize() override;

   protected:
    TableOutput(const std::string& filename, const nlohmann::json& system_config, int interval);
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
    EnergyBalanceOutput(const nlohmann::json& system_config, int interval)
        : TableOutput("energy.txt", system_config, interval) {
    }

   protected:
    const std::vector<std::string>& column_order(const Diagnostics& diag) const override {
        return diag.energyOrder;
    }
};

// === Boundary stats (boundary.txt) ===
class BoundaryOutput : public TableOutput {
   public:
    BoundaryOutput(const nlohmann::json& system_config, int interval)
        : TableOutput("boundary.txt", system_config, interval) {
    }

   protected:
    const std::vector<std::string>& column_order(const Diagnostics& diag) const override {
        return diag.boundaryOrder;
    }
};

// === Console summary (per-step print) ===
class ConsoleOutput : public IDiagnosticOutput {
    const Species& species_;
    double dt_;

   public:
    ConsoleOutput(const Species& species, const nlohmann::json& system_config);
    void output(int timestep, Diagnostics& energy_diag) override;
};

// === 2D field slices (E, B) ===
class Fields2DOutput : public IDiagnosticOutput {
    const Domain& domain_;
    const Field3d& fieldEn_;
    const Field3d& fieldBn_;
    const Field3d& fieldE_external_;
    const Field3d& fieldBInit_;
    std::vector<std::string> field_names_;
    std::vector<double> planes_x_, planes_y_, planes_z_;
    int interval_;

    // Усреднение по окну average_steps_ шагов (0 — выключено)
    int average_steps_ = 0;
    AveragedField2D average_;

    // E/B are sampled as sums at plane nodes, without building a 3D temporary.
    std::pair<const Field3d*, const Field3d*> field_components(const std::string& fname) const;
    void write_averages(int sequence) const;

   public:
    Fields2DOutput(const Domain& domain, const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                   const Field3d& fieldBInit, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
    void finalize() override;
};

// === Base for species 2D outputs (density, current, pressure) ===
// Avoids code duplication for plane parsing and species resolution
class Species2DOutput : public IDiagnosticOutput {
   public:
    Species2DOutput(const Species& species, const Domain& domain, const nlohmann::json& config);
    virtual ~Species2DOutput() = default;
    void output(int timestep, Diagnostics&) override;
    void finalize() override;

   protected:
    const Species& species_;
    const Domain& domain_;
    std::vector<std::string> species_names_;
    std::vector<double> planes_x_, planes_y_, planes_z_;
    int interval_;

    // Усреднение по окну average_steps_ шагов (0 — выключено)
    int average_steps_ = 0;
    AveragedField2D average_;

    const ParticlesArray& require_species(const std::string& name) const;
    void write_averages(int sequence) const;

    // The callback consumes a field immediately, so expensive moments can reuse
    // one temporary buffer for all components and both output modes.
    using FieldVisitor = std::function<void(const Field3d&, const std::string&)>;
    virtual void for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const = 0;
};

// === 2D density slices per species ===
class Density2DOutput : public Species2DOutput {
   public:
    using Species2DOutput::Species2DOutput;

   protected:
    void for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const override;
};

// === 2D current slices per species ===
class Current2DOutput : public Species2DOutput {
   public:
    using Species2DOutput::Species2DOutput;

   protected:
    void for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const override;
};

// === 2D pressure tensor per species (6 components: Prr, Ppp, Pzz, Prp, Prz, Pzp) ===
class Pressure2DOutput : public Species2DOutput {
   public:
    using Species2DOutput::Species2DOutput;

   protected:
    void for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const override;
};

// === Energy spectrum text file per species ===
class EnergySpectrumOutput : public IDiagnosticOutput {
    const Species& species_;
    std::vector<std::string> species_names_;
    int interval_;

   public:
    EnergySpectrumOutput(const Species& species, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
};

// === Checkpoint / recovery ===
class RecoveryOutput : public IDiagnosticOutput {
    const Field3d& fieldEn_;
    const Field3d& fieldBn_;
    const Species& species_;
    int interval_;

   public:
    RecoveryOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Species& species, int interval);
    void output(int timestep, Diagnostics&) override;
};

// === 3D field dump at specific timesteps ===
class Fields3DOutput : public IDiagnosticOutput {
    const Field3d& fieldEn_;
    const Field3d& fieldBn_;
    const Field3d& fieldE_external_;
    const Field3d& fieldBInit_;
    std::vector<std::string> field_names_;
    std::vector<int> timesteps_;

   public:
    Fields3DOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                   const Field3d& fieldBInit, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
};

// === 2D total charge density (sum over all species) ===
class ChargeDensity2DOutput : public IDiagnosticOutput {
    const Species& species_;
    const Domain& domain_;
    std::vector<double> planes_x_, planes_y_, planes_z_;
    int interval_;

   public:
    ChargeDensity2DOutput(const Species& species, const Domain& domain, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
};

// === 2D total current (sum over all charged species) ===
class TotalCurrent2DOutput : public IDiagnosticOutput {
    const Species& species_;
    const Domain& domain_;
    std::vector<double> planes_x_, planes_y_, planes_z_;
    int interval_;

   public:
    TotalCurrent2DOutput(const Species& species, const Domain& domain, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
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
                    const Field3d& fieldBInit, const Species& species, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
};

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
                const Field3d& fieldBInit, const Species& species, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
    void finalize() override;
};

// === Factory: creates IDiagnosticOutput vector from JSON config ===
class OutputFactory {
   public:
    static std::vector<std::unique_ptr<IDiagnosticOutput>> create(const nlohmann::json& diag_config,
                                                                  const Domain& domain, const Species& species,
                                                                  const nlohmann::json& system_config,
                                                                  const Field3d& fieldEn, const Field3d& fieldBn,
                                                                  const Field3d& fieldE_external,
                                                                  const Field3d& fieldBInit);

   private:
    static std::vector<std::unique_ptr<IDiagnosticOutput>> create_from_list(
        const nlohmann::json& outputs, const Domain& domain, const Species& species,
        const nlohmann::json& system_config, const Field3d& fieldEn, const Field3d& fieldBn,
        const Field3d& fieldE_external, const Field3d& fieldBInit);
};

#endif
