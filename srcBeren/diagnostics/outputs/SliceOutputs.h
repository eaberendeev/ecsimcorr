#pragma once

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "AveragedField2D.h"
#include "DiagnosticOutputSettings.h"
#include "IDiagnosticOutput.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"

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
                   const Field3d& fieldBInit, const SliceOutputSettings& settings,
                   const std::vector<std::string>& field_names);
    void output(int timestep, const Diagnostics&) override;
    void finalize() override;
};

// === Base for species 2D outputs (density, current, pressure) ===
// Shares snapshot/average scheduling and reuses one moment buffer per species.
class Species2DOutput : public IDiagnosticOutput {
   public:
    Species2DOutput(const Species& species, const Domain& domain, const SliceOutputSettings& settings,
                    const std::vector<std::string>& species_names);
    virtual ~Species2DOutput() = default;
    void output(int timestep, const Diagnostics&) override;
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

// Raw cylindrical second moment: sum(m*mpw*S(x)*v_a*v_b), including bulk flow.
// Legacy name pressure_2d and file names are preserved. No mean velocity is
// subtracted; this is not the thermal pressure tensor. Components: rr, pp, zz, rp, rz, zp.
class Pressure2DOutput : public Species2DOutput {
   public:
    using Species2DOutput::Species2DOutput;

   protected:
    void for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const override;
};

// === 2D total charge density (sum over all species) ===
class ChargeDensity2DOutput : public IDiagnosticOutput {
    const Species& species_;
    const Domain& domain_;
    std::vector<double> planes_x_, planes_y_, planes_z_;
    int interval_;

   public:
    ChargeDensity2DOutput(const Species& species, const Domain& domain, const DiagnosticPlanes& planes, int interval);
    void output(int timestep, const Diagnostics&) override;
};

// === 2D total current (sum over all charged species) ===
class TotalCurrent2DOutput : public IDiagnosticOutput {
    const Species& species_;
    const Domain& domain_;
    std::vector<double> planes_x_, planes_y_, planes_z_;
    int interval_;

   public:
    TotalCurrent2DOutput(const Species& species, const Domain& domain, const DiagnosticPlanes& planes, int interval);
    void output(int timestep, const Diagnostics&) override;
};
