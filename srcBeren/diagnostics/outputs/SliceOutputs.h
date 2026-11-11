#pragma once

#include <string>
#include <vector>

#include "diagnostics/AveragedField2D.h"
#include "diagnostics/DiagnosticMeasurements.h"
#include "diagnostics/DiagnosticOutputSettings.h"
#include "diagnostics/IDiagnosticOutput.h"
#include "grid/domain.h"
#include "particles/particle_array.h"
#include "utils/containers.h"

// === 2D field slices (E, B) ===
class Fields2DOutput : public IDiagnosticOutput {
    std::vector<OutputField> fields_;
    int interval_;

    // Усреднение по окну average_steps_ шагов (0 — выключено)
    int average_steps_ = 0;
    AveragedField2D average_;

    void write_averages(int sequence) const;

   public:
    Fields2DOutput(const Domain& domain, const SliceOutputSettings& settings, const std::vector<OutputField>& fields);
    void output(int timestep, const Diagnostics&, const DiagnosticMeasurements&) override;
    void finalize() override;
};

// Density is measured on request; current is prepared by the scheme.
// Second moments are exposed as selected planes, including bulk flow.
enum class SpeciesGridQuantity { Density, Current, SecondMoment };

class Species2DOutput : public IDiagnosticOutput {
   public:
    // `species_names` may include guiding-centre species (DiagnosticSources::extra_species).
    Species2DOutput(const Domain& domain, const SliceOutputSettings& settings,
                    const std::vector<std::string>& species_names, SpeciesGridQuantity quantity);
    void request_measurements(int timestep, DiagnosticMeasurements& measured) const override;
    void output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) override;
    void finalize() override;

   private:
    std::vector<std::string> species_names_;
    int interval_;
    int average_steps_;
    SpeciesGridQuantity quantity_;
    AveragedField2D average_;
    void write_averages(int sequence) const;
};

// Sum only charged species; density is scalar and current has three components.
class TotalGrid2DOutput : public IDiagnosticOutput {
    SliceIndices planes_;
    int interval_;
    TotalGridQuantity quantity_;

   public:
    TotalGrid2DOutput(const Domain& domain, const DiagnosticPlanes& planes, int interval, TotalGridQuantity quantity);
    void request_measurements(int timestep, DiagnosticMeasurements& measured) const override;
    void output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) override;
};
