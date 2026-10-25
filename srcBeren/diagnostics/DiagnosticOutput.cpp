#include "DiagnosticOutput.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

#include "ParticlesDiagnostic.h"
#include "Read.h"
#include "interpolation.h"
#include "logger.h"
#include "output_util.h"
#include "recovery.h"
#include "timer.h"
#include "util.h"

// ============================================================
//  Helpers (anonymous namespace)
// ============================================================
namespace {

void parse_planes(const nlohmann::json& config, std::vector<double>& px, std::vector<double>& py,
                  std::vector<double>& pz) {
    if (config.contains("planes")) {
        const auto& p = config["planes"];
        if (p.contains("x"))
            px = p["x"].get<std::vector<double>>();
        if (p.contains("y"))
            py = p["y"].get<std::vector<double>>();
        if (p.contains("z"))
            pz = p["z"].get<std::vector<double>>();
    }
}

int get_interval(const nlohmann::json& config, int default_val = 1) {
    const int interval = config.value("interval", default_val);
    if (interval <= 0) {
        const std::string msg = "Diagnostic output interval must be > 0 (got " + std::to_string(interval) + ")";
        logger::error(msg);
        throw std::runtime_error(msg);
    }
    return interval;
}

// For 2D slices, an absent interval means no instantaneous output.
int get_slice_interval(const nlohmann::json& config) {
    return config.contains("interval") ? get_interval(config) : 0;
}

// Окно усреднения (в шагах); 0 — усреднение выключено, отрицательное — ошибка.
int get_average_steps(const nlohmann::json& config) {
    const int steps = config.value("average_steps", 0);
    if (steps < 0) {
        const std::string msg = "average_steps must be >= 0 (got " + std::to_string(steps) + ")";
        logger::error(msg);
        throw std::runtime_error(msg);
    }
    return steps;
}

// Имя файла (basename) для поля по его короткому имени из конфига.
std::string field_output_basename(const std::string& fname) {
    if (fname == "E")
        return "FieldE";
    if (fname == "B")
        return "FieldB";
    if (fname == "En")
        return "FieldEn";
    if (fname == "Bn")
        return "FieldBn";
    if (fname == "E_ex")
        return "FieldE_external";
    if (fname == "B_init")
        return "FieldBInit";
    return "";
}

std::vector<std::string> resolve_species_names(const Species& all_species, const nlohmann::json& config) {
    std::vector<std::string> names;
    if (!config.contains("species")) {
        for (const auto& kv : all_species) names.push_back(kv.first);
        return names;
    }
    const auto& spec = config["species"];
    if (spec.is_string() && spec.get<std::string>() == "all") {
        for (const auto& kv : all_species) names.push_back(kv.first);
    } else if (spec.is_array()) {
        names = spec.get<std::vector<std::string>>();
    }
    return names;
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
//  AveragedField2D
// ============================================================
void AveragedField2D::set_planes(const std::vector<double>& px, const std::vector<double>& py,
                                 const std::vector<double>& pz, const Domain& domain) {
    const std::array<std::vector<double>, 3> coords{px, py, pz};
    for (int axis = 0; axis < 3; ++axis) {
        positions_[axis].clear();
        for (double coord : coords[axis]) positions_[axis].push_back(int_value(coord / domain.cell_size()[axis]));
    }
}

void AveragedField2D::add(const std::string& key, const Field3d& field, const Field3d* extra) {
    RECORD_TIMER;
    if (field.capacity() == 0)
        return;
    auto [it, inserted] = accum_.try_emplace(key);
    auto& planes = it->second;
    if (inserted) {
        for (int axis = 0; axis < 3; ++axis) {
            for (int position : positions_[axis]) {
                if (position < 0 || position >= field.sizes()[axis])
                    continue;   // Same out-of-range behavior as instantaneous slices.
                auto sizes = field.sizes();
                sizes[axis] = 1;
                planes.push_back({axis, position, Field3d(sizes, field.nd())});
                planes.back().sum.setZero();
            }
        }
    }
    for (auto& plane : planes) {
        const auto sizes = plane.sum.sizes();
        for (int i = 0; i < sizes.x(); ++i)
            for (int j = 0; j < sizes.y(); ++j)
                for (int k = 0; k < sizes.z(); ++k) {
                    Vector3I node(i, j, k);
                    node[plane.axis] = plane.position;
                    for (int d = 0; d < field.nd(); ++d) {
                        double value = field(node.x(), node.y(), node.z(), d);
                        if (extra)
                            value += (*extra)(node.x(), node.y(), node.z(), d);
                        plane.sum(i, j, k, d) += value;
                    }
                }
    }
}

void AveragedField2D::write(const std::string& key, const std::string& basename, int sequence, bool averaged) const {
    RECORD_TIMER;
    const auto it = accum_.find(key);
    if (it == accum_.end() || count == 0)
        return;
    const std::string suffix = std::string(averaged ? "avg" : "") + to_string(sequence, 4);
    for (const auto& plane : it->second) {
        // Preserve the binary layout: three component slots, even for scalar
        // density/pressure files. Unused slots remain zero.
        const auto sizes = plane.sum.sizes();
        const int first = plane.axis == 0 ? 1 : 0;
        const int second = plane.axis == 2 ? 1 : 2;
        const int n1 = sizes[first], n2 = sizes[second];
        std::vector<float> values(3 * n1 * n2, 0.0f);
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j) {
                Vector3I node(0, 0, 0);
                node[first] = i;
                node[second] = j;
                for (int d = 0; d < plane.sum.nd(); ++d)
                    values[output_field_index2d(i, j, d, n1, n2)] =
                        static_cast<float>(plane.sum(node.x(), node.y(), node.z(), d) / count);
            }
        const std::string axis(1, "XYZ"[plane.axis]);
        output_array2d(values, n1, n2, basename + "_Plane" + axis + "_" + to_string(plane.position, 3) + "_" + suffix);
    }
}

size_t AveragedField2D::stored_values() const {
    size_t result = 0;
    for (const auto& [key, planes] : accum_)
        for (const auto& plane : planes) result += plane.sum.capacity();
    return result;
}

std::vector<std::string> AveragedField2D::keys() const {
    std::vector<std::string> result;
    for (const auto& entry : accum_) result.push_back(entry.first);
    return result;
}

long long AveragedField2D::sequence() const {
    return sequence_;
}

void AveragedField2D::reset() {
    RECORD_TIMER;
    for (auto& [key, planes] : accum_)
        for (auto& plane : planes) plane.sum.setZero();
    count = 0;
    ++sequence_;
}

// ============================================================
//  TableOutput  (shared by energy.txt and boundary.txt)
// ============================================================
TableOutput::TableOutput(const std::string& filename, const nlohmann::json& system_config, int interval)
    : filename_(filename), dt_(get_checked<double>(system_config, "Dt")), interval_(interval), header_written_(false) {
    file_.open(filename_);
    if (!file_.is_open()) {
        const std::string msg = "Cannot open " + filename_ + " for writing";
        logger::error(msg);
        throw std::runtime_error(msg);
    }
}

void TableOutput::output(int timestep, Diagnostics& diag) {
    if (timestep % interval_ != 0)
        return;

    std::stringstream ss;
    if (!header_written_) {
        ss << "Time ";
        for (const auto& key : column_order(diag)) {
            ss << key << " ";
        }
        header_written_ = true;
        ss << "\n";
    }
    // 12 digits: energy conservation values are ~1e-13 and need more than the
    // default 6 digits to be useful.
    ss << std::setprecision(12) << timestep * dt_ << " ";
    for (const auto& key : column_order(diag)) {
        ss << diag.energy[key] << " ";
    }
    ss << "\n";
    file_ << ss.str();
    file_.flush();
}

void TableOutput::finalize() {
    if (file_.is_open()) {
        file_.flush();
        file_.close();
    }
}

// ============================================================
//  ConsoleOutput
// ============================================================
ConsoleOutput::ConsoleOutput(const Species& species, const nlohmann::json& system_config)
    : species_(species), dt_(get_checked<double>(system_config, "Dt")) {
}

void ConsoleOutput::output(int timestep, Diagnostics& diag) {
    const double t = timestep * dt_;
    std::ostringstream line;
    line << std::setw(6) << timestep << "  " << std::fixed << std::setprecision(2) << std::setw(10) << t
         << std::scientific << std::setprecision(3) << "  Ef=" << std::setw(10) << diag.energy["energyFieldE"]
         << " Bf=" << std::setw(10) << diag.energy["energyFieldB"];

    for (const auto& kv : species_) {
        const auto& sp = *kv.second;
        const int n = sp.get_total_num_of_particles();
        const auto key = sp.name();
        const double ek = diag.energy.count(key) ? diag.energy.at(key) : 0.0;
        line << "  " << sp.name()[0] << ":" << std::setw(5) << n << " E=" << std::setw(10) << ek;
    }

    line << "  dE=" << std::setw(9) << diag.energy["energyConserve"];
    std::cout << line.str() << "\n";
}

// ============================================================
//  Fields2DOutput
// ============================================================
Fields2DOutput::Fields2DOutput(const Domain& domain, const Field3d& fieldEn, const Field3d& fieldBn,
                               const Field3d& fieldE_external, const Field3d& fieldBInit, const nlohmann::json& config)
    : domain_(domain),
      fieldEn_(fieldEn),
      fieldBn_(fieldBn),
      fieldE_external_(fieldE_external),
      fieldBInit_(fieldBInit),
      interval_(get_slice_interval(config)),
      average_steps_(get_average_steps(config)) {
    if (config.contains("fields")) {
        field_names_ = config["fields"].get<std::vector<std::string>>();
    } else {
        field_names_ = {"E", "B"};
    }
    parse_planes(config, planes_x_, planes_y_, planes_z_);
    average_.set_planes(planes_x_, planes_y_, planes_z_, domain_);
}

std::pair<const Field3d*, const Field3d*> Fields2DOutput::field_components(const std::string& fname) const {
    if (fname == "E")
        return {&fieldEn_, &fieldE_external_};
    if (fname == "B")
        return {&fieldBn_, &fieldBInit_};
    if (fname == "En")
        return {&fieldEn_, nullptr};
    if (fname == "Bn")
        return {&fieldBn_, nullptr};
    if (fname == "E_ex")
        return {&fieldE_external_, nullptr};
    if (fname == "B_init")
        return {&fieldBInit_, nullptr};
    return {nullptr, nullptr};
}

void Fields2DOutput::write_averages(int sequence) const {
    const std::string base = ".//Fields//Diag2D//";
    for (const auto& fname : field_names_) {
        const std::string name = field_output_basename(fname);
        if (name.empty())
            continue;
        average_.write(fname, base + name, sequence);
    }
}

void Fields2DOutput::output(int timestep, Diagnostics&) {
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
    for (const auto& fname : field_names_) {
        const auto [field, extra] = field_components(fname);
        if (!field)
            continue;
        if (average_steps_ > 0)
            average_.add(fname, *field, extra);
        if (write_now) {
            instantaneous.add(fname, *field, extra);
            instantaneous.write(fname, base + field_output_basename(fname), timestep / interval_, false);
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
Species2DOutput::Species2DOutput(const Species& species, const Domain& domain, const nlohmann::json& config)
    : species_(species),
      domain_(domain),
      interval_(get_slice_interval(config)),
      average_steps_(get_average_steps(config)) {
    species_names_ = resolve_species_names(species, config);
    parse_planes(config, planes_x_, planes_y_, planes_z_);
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

void Species2DOutput::output(int timestep, Diagnostics&) {
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

// ============================================================
//  EnergySpectrumOutput
// ============================================================
EnergySpectrumOutput::EnergySpectrumOutput(const Species& species, const nlohmann::json& config)
    : species_(species), interval_(get_interval(config)) {
    species_names_ = resolve_species_names(species, config);
}

void EnergySpectrumOutput::output(int timestep, Diagnostics&) {
    if (timestep % interval_ != 0)
        return;

    for (const auto& name : species_names_) {
        const auto it = species_.find(name);
        if (it == species_.end()) {
            const std::string msg = "EnergySpectrumOutput: species '" + name + "' not found";
            logger::error(msg);
            throw std::runtime_error(msg);
        }

        const auto& sp = *it->second;
        EnergySpectrum spectrum = sp.calculate_energy_spectrum();

        std::ostringstream filename;
        filename << ".//Particles//" << name << "//energy_spectrum_" << std::setw(6) << std::setfill('0') << timestep
                 << ".txt";

        std::ofstream out(filename.str());
        if (!out.is_open()) {
            const std::string msg = "Cannot open " + filename.str() + " for writing";
            logger::error(msg);
            throw std::runtime_error(msg);
        }

        const int nbins = spectrum.spectrum.size();
        const double bin_width = (nbins > 0) ? (spectrum.maxEnergy - spectrum.minEnergy) / nbins : 0.0;
        for (int i = 0; i < nbins; ++i) {
            const double ec = spectrum.minEnergy + (i + 0.5) * bin_width;
            out << ec << " " << spectrum.spectrum[i] << "\n";
        }
    }
}

// ============================================================
//  RecoveryOutput
// ============================================================
RecoveryOutput::RecoveryOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Species& species, int interval)
    : fieldEn_(fieldEn), fieldBn_(fieldBn), species_(species), interval_(interval) {
    create_directory(".//Recovery");
    create_directory(".//Recovery//Fields");
    create_directory(".//Recovery//Particles");
    for (const auto& kv : species_) {
        create_directory(".//Recovery//Particles//" + kv.second->name());
    }
}

void RecoveryOutput::output(int timestep, Diagnostics&) {
    if (timestep % interval_ != 0)
        return;
    for (const auto& kv : species_) {
        write_particles_to_recovery(kv.second.get(), timestep, interval_);
    }
    write_fields_to_recovery(fieldEn_, fieldBn_, timestep, interval_);
}

// ============================================================
//  Fields3DOutput
// ============================================================
Fields3DOutput::Fields3DOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                               const Field3d& fieldBInit, const nlohmann::json& config)
    : fieldEn_(fieldEn), fieldBn_(fieldBn), fieldE_external_(fieldE_external), fieldBInit_(fieldBInit) {
    if (config.contains("fields")) {
        field_names_ = config["fields"].get<std::vector<std::string>>();
    } else {
        field_names_ = {"E", "B"};
    }
    if (config.contains("timesteps")) {
        timesteps_ = config["timesteps"].get<std::vector<int>>();
    }
    create_directory(".//Fields//Diag3D");
}

void Fields3DOutput::output(int timestep, Diagnostics&) {
    if (std::find(timesteps_.begin(), timesteps_.end(), timestep) == timesteps_.end())
        return;

    std::string dir = ".//Fields//Diag3D//";
    auto write_field = [&](const Field3d& field, const std::string& name) {
        std::string fname = dir + name + "_" + std::to_string(timestep) + ".3d";
        write_field_to_file(fname, field);
    };

    for (const auto& fname : field_names_) {
        if (fname == "E") {
            Field3d fieldEFull = fieldEn_ + fieldE_external_;
            write_field(fieldEFull, "FieldE");
        } else if (fname == "B") {
            Field3d fieldBFull(fieldBn_.sizes(), fieldBn_.nd());
            fieldBFull.data() = fieldBn_.data() + fieldBInit_.data();
            write_field(fieldBFull, "FieldB");
        } else if (fname == "En") {
            write_field(fieldEn_, "FieldEn");
        } else if (fname == "Bn") {
            write_field(fieldBn_, "FieldBn");
        } else if (fname == "E_ex") {
            write_field(fieldE_external_, "FieldE_external");
        }
    }
}

// ============================================================
//  ChargeDensity2DOutput
// ============================================================
ChargeDensity2DOutput::ChargeDensity2DOutput(const Species& species, const Domain& domain, const nlohmann::json& config)
    : species_(species), domain_(domain), interval_(get_interval(config)) {
    parse_planes(config, planes_x_, planes_y_, planes_z_);
}

void ChargeDensity2DOutput::output(int timestep, Diagnostics&) {
    if (timestep % interval_ != 0)
        return;
    int seq = std::max(0, timestep / interval_);

    Field3d totalDensity(domain_.size(), 1);
    totalDensity.setZero();
    for (const auto& kv : species_) {
        if (!kv.second->is_neutral())
            totalDensity.data() += kv.second->densityOnGrid.data();
    }

    write_slices(totalDensity, ".//Fields//Diag2D//ChargeDensity", planes_x_, planes_y_, planes_z_, domain_, seq);
}

// ============================================================
//  TotalCurrent2DOutput
// ============================================================
TotalCurrent2DOutput::TotalCurrent2DOutput(const Species& species, const Domain& domain, const nlohmann::json& config)
    : species_(species), domain_(domain), interval_(get_interval(config)) {
    parse_planes(config, planes_x_, planes_y_, planes_z_);
}

void TotalCurrent2DOutput::output(int timestep, Diagnostics&) {
    if (timestep % interval_ != 0)
        return;
    int seq = std::max(0, timestep / interval_);

    Field3d totalCurrent(domain_.size(), 3);
    totalCurrent.setZero();
    for (const auto& kv : species_) {
        if (!kv.second->is_neutral())
            totalCurrent.data() += kv.second->currentOnGrid.data();
    }

    write_slices(totalCurrent, ".//Fields//Diag2D//TotalCurrent", planes_x_, planes_y_, planes_z_, domain_, seq);
}

// ============================================================
//  Full3DAllOutput  — full dump at given timesteps
// ============================================================
Full3DAllOutput::Full3DAllOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                                 const Field3d& fieldBInit, const Species& species, const nlohmann::json& config)
    : fieldEn_(fieldEn),
      fieldBn_(fieldBn),
      fieldE_external_(fieldE_external),
      fieldBInit_(fieldBInit),
      species_(species) {
    if (config.contains("timesteps")) {
        timesteps_ = config["timesteps"].get<std::vector<int>>();
    }
    create_directory(".//Full3D");
    create_directory(".//Full3D//Fields");
    create_directory(".//Full3D//Particles");
    for (const auto& kv : species_) {
        create_directory(".//Full3D//Particles//" + kv.second->name());
    }
}

void Full3DAllOutput::output(int timestep, Diagnostics&) {
    if (std::find(timesteps_.begin(), timesteps_.end(), timestep) == timesteps_.end())
        return;

    std::string tstr = std::to_string(timestep);

    auto dump = [&](const Field3d& f, const std::string& path) { write_field_to_file(path, f); };

    dump(fieldEn_, ".//Full3D//Fields//En_" + tstr + ".3d");
    dump(fieldBn_, ".//Full3D//Fields//Bn_" + tstr + ".3d");
    dump(fieldE_external_, ".//Full3D//Fields//E_ex_" + tstr + ".3d");
    dump(fieldBInit_, ".//Full3D//Fields//B_init_" + tstr + ".3d");

    Field3d eFull(fieldEn_.sizes(), fieldEn_.nd());
    eFull.data() = fieldEn_.data() + fieldE_external_.data();
    dump(eFull, ".//Full3D//Fields//E_" + tstr + ".3d");

    Field3d bFull(fieldBn_.sizes(), fieldBn_.nd());
    bFull.data() = fieldBn_.data() + fieldBInit_.data();
    dump(bFull, ".//Full3D//Fields//B_" + tstr + ".3d");

    for (const auto& kv : species_) {
        const auto& sp = *kv.second;
        std::string dir = ".//Full3D//Particles//" + sp.name() + "//";
        dump(sp.densityOnGrid, dir + "Density_" + tstr + ".3d");
        dump(sp.currentOnGrid, dir + "Current_" + tstr + ".3d");
    }
}

// ============================================================
//  ProbeOutput  — interpolate fields at probe points, write probes.txt
// ============================================================
namespace {

double interpolate_scalar_cell(const Field3d& field, int ix, int iy, int iz, double wx, double wy, double wz) {
    const double w000 = (1 - wx) * (1 - wy) * (1 - wz);
    const double w001 = (1 - wx) * (1 - wy) * wz;
    const double w010 = (1 - wx) * wy * (1 - wz);
    const double w011 = (1 - wx) * wy * wz;
    const double w100 = wx * (1 - wy) * (1 - wz);
    const double w101 = wx * (1 - wy) * wz;
    const double w110 = wx * wy * (1 - wz);
    const double w111 = wx * wy * wz;
    const Field3d& f = field;
    return w000 * f(ix, iy, iz, 0) + w001 * f(ix, iy, iz + 1, 0) + w010 * f(ix, iy + 1, iz, 0) +
           w011 * f(ix, iy + 1, iz + 1, 0) + w100 * f(ix + 1, iy, iz, 0) + w101 * f(ix + 1, iy, iz + 1, 0) +
           w110 * f(ix + 1, iy + 1, iz, 0) + w111 * f(ix + 1, iy + 1, iz + 1, 0);
}

Vector3R interpolate_vector_at(const Field3d& field, const Vector3R& world, const Domain& domain, bool is_b) {
    const Vector3R cell = domain.cell_size();
    if (is_b) {
        return interpolateB_linear(field, Vector3R(world.x() / cell.x(), world.y() / cell.y(), world.z() / cell.z()));
    } else {
        return interpolateE_linear(field, Vector3R(world.x() / cell.x(), world.y() / cell.y(), world.z() / cell.z()));
    }
}

double interpolate_scalar_at(const Field3d& field, const Vector3R& world, const Domain& domain) {
    const Vector3R cell = domain.cell_size();
    const double x = world.x() / cell.x() + GHOST_CELLS;
    const double y = world.y() / cell.y() + GHOST_CELLS;
    const double z = world.z() / cell.z() + GHOST_CELLS;
    const int ix = int(x), iy = int(y), iz = int(z);
    return interpolate_scalar_cell(field, ix, iy, iz, x - ix, y - iy, z - iz);
}

}   // namespace

ProbeOutput::ProbeOutput(const Domain& domain, const Field3d& fieldEn, const Field3d& fieldBn,
                         const Field3d& fieldE_external, const Field3d& fieldBInit, const Species& species,
                         const nlohmann::json& config)
    : domain_(domain),
      fieldEn_(fieldEn),
      fieldBn_(fieldBn),
      fieldE_external_(fieldE_external),
      fieldBInit_(fieldBInit),
      species_(species),
      interval_(get_interval(config)),
      header_written_(false) {
    species_names_ = resolve_species_names(species, config);
    if (config.contains("points")) {
        for (const auto& pt : config["points"]) {
            if (pt.is_array() && pt.size() == 3) {
                points_.push_back(Vector3R(pt[0].get<double>(), pt[1].get<double>(), pt[2].get<double>()));
            }
        }
    }
    create_directory(".//Probes");
    fProbes_.open(".//Probes//probes.txt");
}

void ProbeOutput::output(int timestep, Diagnostics&) {
    if (timestep % interval_ != 0)
        return;
    if (points_.empty())
        return;

    if (!header_written_) {
        fProbes_ << "# Time x y z";
        fProbes_ << " Ex Ey Ez Bx By Bz";
        fProbes_ << " Enx Eny Enz Bnx Bny Bnz";
        fProbes_ << " E_ex_x E_ex_y E_ex_z B_init_x B_init_y B_init_z";
        for (const auto& name : species_names_) fProbes_ << " dens_" << name;
        fProbes_ << "\n";
        header_written_ = true;
    }

    auto field_e_full = [&]() {
        Field3d f(fieldEn_.sizes(), fieldEn_.nd());
        f.data() = fieldEn_.data() + fieldE_external_.data();
        return f;
    };
    auto field_b_full = [&]() {
        Field3d f(fieldBn_.sizes(), fieldBn_.nd());
        f.data() = fieldBn_.data() + fieldBInit_.data();
        return f;
    };

    for (const auto& pt : points_) {
        fProbes_ << timestep << " " << pt.x() << " " << pt.y() << " " << pt.z();

        const Field3d eTmp = field_e_full();
        const Field3d bTmp = field_b_full();
        fProbes_ << std::scientific << std::setprecision(6);
        Vector3R w = interpolate_vector_at(eTmp, pt, domain_, false);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(bTmp, pt, domain_, true);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(fieldEn_, pt, domain_, false);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        {
            Field3d bOnly(fieldBn_.sizes(), fieldBn_.nd());
            bOnly.data() = fieldBn_.data();
            w = interpolate_vector_at(bOnly, pt, domain_, true);
        }
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(fieldE_external_, pt, domain_, false);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();
        w = interpolate_vector_at(fieldBInit_, pt, domain_, true);
        fProbes_ << " " << w.x() << " " << w.y() << " " << w.z();

        for (const auto& name : species_names_) {
            const auto it = species_.find(name);
            double d = 0;
            if (it != species_.end()) {
                d = interpolate_scalar_at(it->second->densityOnGrid, pt, domain_);
            }
            fProbes_ << " " << d;
        }
        fProbes_ << "\n";
    }
    fProbes_.flush();
}

void ProbeOutput::finalize() {
    if (fProbes_.is_open()) {
        fProbes_.flush();
        fProbes_.close();
    }
}

// ============================================================
//  OutputFactory
// ============================================================
std::vector<std::unique_ptr<IDiagnosticOutput>> OutputFactory::create(const nlohmann::json& diag_config,
                                                                      const Domain& domain, const Species& species,
                                                                      const nlohmann::json& system_config,
                                                                      const Field3d& fieldEn, const Field3d& fieldBn,
                                                                      const Field3d& fieldE_external,
                                                                      const Field3d& fieldBInit) {
    if (!diag_config.contains("outputs"))
        return {};
    if (!diag_config["outputs"].is_array()) {
        throw std::runtime_error(
            "OutputFactory: missing 'outputs' array in diagnostics config. "
            "Update your gen_config.py to use the new format.");
    }
    return create_from_list(diag_config["outputs"], domain, species, system_config, fieldEn, fieldBn, fieldE_external,
                            fieldBInit);
}

std::vector<std::unique_ptr<IDiagnosticOutput>> OutputFactory::create_from_list(
    const nlohmann::json& outputs, const Domain& domain, const Species& species, const nlohmann::json& system_config,
    const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external, const Field3d& fieldBInit) {
    std::vector<std::unique_ptr<IDiagnosticOutput>> result;

    for (const auto& cfg : outputs) {
        std::string type = cfg.value("type", "");

        if (type == "energy_balance") {
            result.push_back(std::make_unique<EnergyBalanceOutput>(system_config, get_interval(cfg)));
        } else if (type == "boundary_stats") {
            result.push_back(std::make_unique<BoundaryOutput>(system_config, get_interval(cfg)));
        } else if (type == "console_summary") {
            result.push_back(std::make_unique<ConsoleOutput>(species, system_config));
        } else if (type == "fields_2d") {
            result.push_back(
                std::make_unique<Fields2DOutput>(domain, fieldEn, fieldBn, fieldE_external, fieldBInit, cfg));
        } else if (type == "density_2d") {
            result.push_back(std::make_unique<Density2DOutput>(species, domain, cfg));
        } else if (type == "current_2d") {
            result.push_back(std::make_unique<Current2DOutput>(species, domain, cfg));
        } else if (type == "pressure_2d") {
            result.push_back(std::make_unique<Pressure2DOutput>(species, domain, cfg));
        } else if (type == "energy_spectrum") {
            result.push_back(std::make_unique<EnergySpectrumOutput>(species, cfg));
        } else if (type == "recovery") {
            // Recovery is optional: absent interval disables it, but an
            // explicitly configured interval must be positive.
            if (cfg.contains("interval")) {
                const int interval = get_interval(cfg);
                result.push_back(std::make_unique<RecoveryOutput>(fieldEn, fieldBn, species, interval));
            }
        } else if (type == "fields_3d") {
            result.push_back(std::make_unique<Fields3DOutput>(fieldEn, fieldBn, fieldE_external, fieldBInit, cfg));
        } else if (type == "charge_density_2d") {
            result.push_back(std::make_unique<ChargeDensity2DOutput>(species, domain, cfg));
        } else if (type == "total_current_2d") {
            result.push_back(std::make_unique<TotalCurrent2DOutput>(species, domain, cfg));
        } else if (type == "full3d_all") {
            result.push_back(
                std::make_unique<Full3DAllOutput>(fieldEn, fieldBn, fieldE_external, fieldBInit, species, cfg));
        } else if (type == "probes") {
            result.push_back(
                std::make_unique<ProbeOutput>(domain, fieldEn, fieldBn, fieldE_external, fieldBInit, species, cfg));
        } else {
            const std::string msg = "OutputFactory: unknown output type '" + type + "'";
            logger::error(msg);
            throw std::runtime_error(msg);
        }
    }
    return result;
}
