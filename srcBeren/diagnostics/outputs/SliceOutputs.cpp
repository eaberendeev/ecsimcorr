#include "SliceOutputs.h"

#include <algorithm>
#include <stdexcept>

#include "ParticlesDiagnostic.h"
#include "SliceFiles.h"
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

void write_slices(const Field3d& field, const std::string& basename, const std::vector<double>& px,
                  const std::vector<double>& py, const std::vector<double>& pz, const Domain& domain, int seqnum,
                  bool averaged = false) {
    const Vector3R cell = domain.cell_size();
    // усреднённые файлы: <basename>_Plane<axis>_<pos>_avg<NNNN>
    const std::string sNumber = std::string(averaged ? "avg" : "") + to_string(seqnum, 4);
    const Vector3I start = {0, 0, 0};
    const Vector3I end = field.sizes();

    auto write_axis = [&](const std::vector<double>& coords, int axis, double dr) {
        for (double coord : coords) {
            const int pos = int_value(coord / dr);
            output_field_plane(field, start, end, pos, axis, field.nd(), basename, sNumber);
        }
    };
    write_axis(px, 0, cell.x());
    write_axis(py, 1, cell.y());
    write_axis(pz, 2, cell.z());
}

}   // namespace

// ============================================================
//  Fields2DOutput
// ============================================================
Fields2DOutput::Fields2DOutput(const Domain& domain, const SliceOutputSettings& settings,
                               const std::vector<OutputField>& fields)
    : domain_(domain),
      fields_(fields),
      planes_x_(settings.planes.x),
      planes_y_(settings.planes.y),
      planes_z_(settings.planes.z),
      interval_(settings.interval),
      average_steps_(settings.average_steps) {
    prepare_field_directory("Diag2D");
    average_.set_planes(planes_x_, planes_y_, planes_z_, domain_);
}

void Fields2DOutput::write_averages(int sequence) const {
    const std::string base = ".//Fields//Diag2D//";
    for (const auto& selected : fields_) average_.write(selected.basename, base + selected.basename, sequence);
}

void Fields2DOutput::output(int timestep, const Diagnostics&) {
    RECORD_TIMER;
    const std::string base = ".//Fields//Diag2D//";

    const bool write_now = interval_ > 0 && timestep % interval_ == 0;
    if (!write_now && average_steps_ == 0)
        return;
    AveragedField2D instantaneous;
    if (write_now) {
        instantaneous.set_planes(planes_x_, planes_y_, planes_z_, domain_);
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
//  Species2DOutput  (base for density, current, pressure)
// ============================================================
Species2DOutput::Species2DOutput(const Species& species, const Domain& domain, const SliceOutputSettings& settings,
                                 const std::vector<std::string>& species_names)
    : species_(species),
      domain_(domain),
      species_names_(species_names),
      planes_x_(settings.planes.x),
      planes_y_(settings.planes.y),
      planes_z_(settings.planes.z),
      interval_(settings.interval),
      average_steps_(settings.average_steps) {
    prepare_species_directories(species_names_);
    average_.set_planes(planes_x_, planes_y_, planes_z_, domain_);
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

void Species2DOutput::output(int timestep, const Diagnostics&) {
    RECORD_TIMER;
    const bool write_now = interval_ > 0 && timestep % interval_ == 0;
    if (!write_now && average_steps_ == 0)
        return;
    for (const auto& name : species_names_) {
        for_each_field(require_species(name), [&](const Field3d& field, const std::string& basename) {
            if (average_steps_ > 0)
                average_.add(basename, field);
            if (write_now)
                write_slices(field, basename, planes_x_, planes_y_, planes_z_, domain_, timestep / interval_);
        });
    }
    if (average_steps_ > 0) {
        ++average_.count;
        if (average_.count >= average_steps_) {
            write_averages(static_cast<int>(average_.sequence()));
            average_.reset();
        }
    }
}

void Species2DOutput::finalize() {
    // Дописываем неполное окно усреднения в конце прогона.
    if (average_.count > 0) {
        write_averages(static_cast<int>(average_.sequence()));
        average_.reset();
    }
}

// ============================================================
//  Density2DOutput
// ============================================================
void Density2DOutput::for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const {
    visit(sp.densityOnGrid, ".//Particles//" + sp.name() + "//Diag2D//Density");
}

void Current2DOutput::for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const {
    visit(sp.currentOnGrid, ".//Particles//" + sp.name() + "//Diag2D//Current");
}

void Pressure2DOutput::for_each_field(const ParticlesArray& sp, const FieldVisitor& visit) const {
    RECORD_TIMER;
    if (sp.is_neutral())
        return;
    const std::string base = ".//Particles//" + sp.name() + "//Diag2D//";
    Field3d moment(domain_.size(), 1);
    const auto component = [&](const char* name, auto v1, auto v2) {
        sp.calculate_pressure_component(moment, v1, v2);
        visit(moment, base + name);
    };
    component("Prr", radial_velocity, radial_velocity);
    component("Ppp", phi_velocity, phi_velocity);
    component("Pzz", z_velocity, z_velocity);
    component("Prp", radial_velocity, phi_velocity);
    component("Prz", radial_velocity, z_velocity);
    component("Pzp", z_velocity, phi_velocity);
}

TotalGrid2DOutput::TotalGrid2DOutput(const Species& species, const Domain& domain, const DiagnosticPlanes& planes,
                                     int interval, TotalGridQuantity quantity)
    : species_(species),
      domain_(domain),
      planes_x_(planes.x),
      planes_y_(planes.y),
      planes_z_(planes.z),
      interval_(interval),
      quantity_(quantity) {
    prepare_field_directory("Diag2D");
}

void TotalGrid2DOutput::output(int timestep, const Diagnostics&) {
    if (timestep % interval_ != 0)
        return;
    const bool charge = quantity_ == TotalGridQuantity::ChargeDensity;
    Field3d total(domain_.size(), charge ? 1 : 3);
    total.setZero();
    for (const auto& kv : species_) {
        const auto& sp = *kv.second;
        if (sp.is_neutral())
            continue;
        const auto& field = charge ? sp.densityOnGrid : sp.currentOnGrid;
        total.data() += field.data();
    }
    const std::string basename = charge ? "ChargeDensity" : "TotalCurrent";
    write_slices(total, ".//Fields//Diag2D//" + basename, planes_x_, planes_y_, planes_z_, domain_,
                 std::max(0, timestep / interval_));
}
