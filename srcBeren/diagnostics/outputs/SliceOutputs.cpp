#include "SliceOutputs.h"

#include <algorithm>
#include <stdexcept>

#include "logger.h"
#include "timer.h"
#include "util.h"

// ============================================================
//  Helpers (anonymous namespace)
// ============================================================
namespace {

void prepare_field_directory(const std::string& subdirectory) {
    ensure_directory(".//Fields");
    ensure_directory(".//Fields//" + subdirectory);
}

void prepare_species_directories(const std::vector<std::string>& names) {
    ensure_directory(".//Particles");
    for (const auto& name : names) {
        ensure_directory(".//Particles//" + name);
        ensure_directory(".//Particles//" + name + "//Diag2D");
    }
}

}   // namespace

// ============================================================
//  Fields2DOutput
// ============================================================
Fields2DOutput::Fields2DOutput(const Domain& domain, const SliceOutputSettings& settings,
                               const std::vector<OutputField>& fields)
    : fields_(fields), interval_(settings.interval), average_steps_(settings.average_steps) {
    prepare_field_directory("Diag2D");
    average_.set_planes(settings.planes.x, settings.planes.y, settings.planes.z, domain);
}

void Fields2DOutput::write_averages(int sequence) const {
    const std::string base = ".//Fields//Diag2D//";
    for (const auto& selected : fields_) average_.write(selected.basename, base + selected.basename, sequence);
}

void Fields2DOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements&) {
    RECORD_TIMER;
    const std::string base = ".//Fields//Diag2D//";

    const bool write_now = interval_ > 0 && timestep % interval_ == 0;
    if (!write_now && average_steps_ == 0)
        return;
    AveragedField2D instantaneous;
    if (write_now) {
        instantaneous.set_planes(average_.planes());
        instantaneous.count = 1;
    }
    for (const auto& selected : fields_) {
        if (average_steps_ > 0)
            average_.add(selected.basename, *selected.field, selected.extra);
        if (write_now) {
            instantaneous.add(selected.basename, *selected.field, selected.extra);
            instantaneous.write(selected.basename, base + selected.basename, timestep / interval_, false);
        }
    }
    if (average_steps_ > 0)
        ++average_.count;

    // Периодически пишем среднее и обнуляем аккумулятор.
    if (average_steps_ > 0 && average_.count >= average_steps_) {
        write_averages(static_cast<int>(average_.sequence()));
        average_.reset();
    }
}

void Fields2DOutput::finalize() {
    // Дописываем неполное окно усреднения в конце прогона.
    if (average_.count > 0) {
        write_averages(static_cast<int>(average_.sequence()));
        average_.reset();
    }
}

// ============================================================
//  Species2DOutput  (density, current, raw second moment)
// ============================================================
Species2DOutput::Species2DOutput(const Species& species, const Domain& domain, const SliceOutputSettings& settings,
                                 const std::vector<std::string>& species_names, SpeciesGridQuantity quantity)
    : species_(species),
      species_names_(species_names),
      interval_(settings.interval),
      average_steps_(settings.average_steps),
      quantity_(quantity) {
    prepare_species_directories(species_names_);
    average_.set_planes(settings.planes.x, settings.planes.y, settings.planes.z, domain);
}

void Species2DOutput::request_measurements(int timestep, DiagnosticMeasurements& measured) const {
    if (quantity_ == SpeciesGridQuantity::Current ||
        (average_steps_ == 0 && (interval_ == 0 || timestep % interval_ != 0)))
        return;
    for (const auto& name : species_names_)
        if (quantity_ == SpeciesGridQuantity::Density)
            measured.request_density(name);
        else
            measured.request_second_moment(name, average_.planes());
}

const ParticlesArray& Species2DOutput::require_species(const std::string& name) const {
    const auto it = species_.find(name);
    if (it == species_.end()) {
        const std::string msg = "Species2DOutput: species '" + name + "' not found";
        logger::error(msg);
        throw std::runtime_error(msg);
    }
    return *it->second;
}

void Species2DOutput::write_averages(int sequence) const {
    for (const auto& basename : average_.keys()) {
        average_.write(basename, basename, sequence);
    }
}

void Species2DOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) {
    RECORD_TIMER;
    const bool write_now = interval_ > 0 && timestep % interval_ == 0;
    if (!write_now && average_steps_ == 0)
        return;
    AveragedField2D instantaneous;
    if (write_now) {
        instantaneous.set_planes(average_.planes());
        instantaneous.count = 1;
    }
    const auto sample = [&](const auto& field, const std::string& basename) {
        if (average_steps_ > 0)
            average_.add(basename, field);
        if (write_now) {
            instantaneous.add(basename, field);
            instantaneous.write(basename, basename, timestep / interval_, false);
        }
    };
    for (const auto& name : species_names_) {
        const std::string base = ".//Particles//" + name + "//Diag2D//";
        if (quantity_ == SpeciesGridQuantity::SecondMoment) {
            const auto& components = measured.second_moment(name);
            for (size_t component = 0; component < components.size(); ++component)
                if (!components[component].empty())
                    sample(components[component], base + second_moment_names[component]);
        } else {
            const auto& species = require_species(name);
            const bool density = quantity_ == SpeciesGridQuantity::Density;
            sample(density ? measured.density(name) : species.currentOnGrid, base + (density ? "Density" : "Current"));
        }
    }
    if (average_steps_ > 0 && ++average_.count >= average_steps_) {
        write_averages(static_cast<int>(average_.sequence()));
        average_.reset();
    }
}

void Species2DOutput::finalize() {
    // Дописываем неполное окно усреднения в конце прогона.
    if (average_.count > 0) {
        write_averages(static_cast<int>(average_.sequence()));
        average_.reset();
    }
}

TotalGrid2DOutput::TotalGrid2DOutput(const Domain& domain, const DiagnosticPlanes& planes, int interval,
                                     TotalGridQuantity quantity)
    : planes_(slice_indices(planes.x, planes.y, planes.z, domain)), interval_(interval), quantity_(quantity) {
    prepare_field_directory("Diag2D");
}

void TotalGrid2DOutput::request_measurements(int timestep, DiagnosticMeasurements& measured) const {
    if (timestep % interval_ == 0)
        measured.request_total(quantity_, planes_);
}

void TotalGrid2DOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) {
    if (timestep % interval_ != 0)
        return;
    AveragedField2D instantaneous;
    instantaneous.set_planes(planes_);
    instantaneous.count = 1;
    const std::string basename = quantity_ == TotalGridQuantity::ChargeDensity ? "ChargeDensity" : "TotalCurrent";
    instantaneous.add(basename, measured.total(quantity_));
    instantaneous.write(basename, ".//Fields//Diag2D//" + basename, std::max(0, timestep / interval_), false);
}
