#pragma once
#ifndef DIAGNOSTIC_OUTPUT_H_
#define DIAGNOSTIC_OUTPUT_H_

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
};

// === Energy balance (energy.txt) ===
class EnergyBalanceOutput : public IDiagnosticOutput {
    std::ofstream fEnergy_;
    double dt_;
    int interval_;
    bool header_written_;

   public:
    EnergyBalanceOutput(const nlohmann::json& system_config, int interval);
    void output(int timestep, Diagnostics& energy_diag) override;
};

// === Boundary stats (boundary.txt) ===
class BoundaryOutput : public IDiagnosticOutput {
    std::ofstream fBoundary_;
    double dt_;
    int interval_;
    bool header_written_;

   public:
    BoundaryOutput(const nlohmann::json& system_config, int interval);
    void output(int timestep, Diagnostics& energy_diag) override;
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

   public:
    Fields2DOutput(const Domain& domain, const Field3d& fieldEn, const Field3d& fieldBn,
                   const Field3d& fieldE_external, const Field3d& fieldBInit,
                   const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
};

// === Base for species 2D outputs (density, current, pressure) ===
// Avoids code duplication for plane parsing and species resolution
class Species2DOutput : public IDiagnosticOutput {
   public:
    Species2DOutput(const Species& species, const Domain& domain, const nlohmann::json& config);
    virtual ~Species2DOutput() = default;
    void output(int timestep, Diagnostics&) override;

   protected:
    const Species& species_;
    const Domain& domain_;
    std::vector<std::string> species_names_;
    std::vector<double> planes_x_, planes_y_, planes_z_;
    int interval_;

    // Subclasses implement this to write what they need for one species
    virtual void write_species(const ParticlesArray& sp, int timestep, int seqnum) const = 0;
};

// === 2D density slices per species ===
class Density2DOutput : public Species2DOutput {
   public:
    using Species2DOutput::Species2DOutput;

   protected:
    void write_species(const ParticlesArray& sp, int timestep, int seqnum) const override;
};

// === 2D current slices per species ===
class Current2DOutput : public Species2DOutput {
   public:
    using Species2DOutput::Species2DOutput;

   protected:
    void write_species(const ParticlesArray& sp, int timestep, int seqnum) const override;
};

// === 2D pressure tensor per species (6 components: Prr, Ppp, Pzz, Prp, Prz, Pzp) ===
class Pressure2DOutput : public Species2DOutput {
   public:
    using Species2DOutput::Species2DOutput;

   protected:
    void write_species(const ParticlesArray& sp, int timestep, int seqnum) const override;
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
    Fields3DOutput(const Field3d& fieldEn, const Field3d& fieldBn,
                   const Field3d& fieldE_external, const Field3d& fieldBInit,
                   const nlohmann::json& config);
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
    Full3DAllOutput(const Field3d& fieldEn, const Field3d& fieldBn,
                    const Field3d& fieldE_external, const Field3d& fieldBInit,
                    const Species& species, const nlohmann::json& config);
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
    ProbeOutput(const Domain& domain, const Field3d& fieldEn, const Field3d& fieldBn,
                const Field3d& fieldE_external, const Field3d& fieldBInit,
                const Species& species, const nlohmann::json& config);
    void output(int timestep, Diagnostics&) override;
};

// === Factory: creates IDiagnosticOutput vector from JSON config ===
class OutputFactory {
   public:
    static std::vector<std::unique_ptr<IDiagnosticOutput>> create(
        const nlohmann::json& diag_config, const Domain& domain, const Species& species,
        const nlohmann::json& system_config, const Field3d& fieldEn, const Field3d& fieldBn,
        const Field3d& fieldE_external, const Field3d& fieldBInit);

   private:
    static std::vector<std::unique_ptr<IDiagnosticOutput>> create_from_list(
        const nlohmann::json& outputs, const Domain& domain, const Species& species,
        const nlohmann::json& system_config, const Field3d& fieldEn, const Field3d& fieldBn,
        const Field3d& fieldE_external, const Field3d& fieldBInit);
};

#endif
