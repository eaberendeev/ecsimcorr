#include "DiagnosticMeasurements.h"

#include <stdexcept>

#include "IDiagnosticOutput.h"
#include "ParticlesDiagnostic.h"
#include "interpolation.h"

namespace {
bool has_planes(const SliceIndices& planes) {
    return !planes[0].empty() || !planes[1].empty() || !planes[2].empty();
}

double sample_density(const Field3d& field, const Vector3R& point, const Domain& domain) {
    const auto cell = domain.to_cell_coordinates(point);
    const double x = cell.x() + GHOST_CELLS, y = cell.y() + GHOST_CELLS, z = cell.z() + GHOST_CELLS;
    const int i = int(x), j = int(y), k = int(z);
    const double wx = x - i, wy = y - j, wz = z - k;
    const double w000 = (1 - wx) * (1 - wy) * (1 - wz), w001 = (1 - wx) * (1 - wy) * wz;
    const double w010 = (1 - wx) * wy * (1 - wz), w011 = (1 - wx) * wy * wz;
    const double w100 = wx * (1 - wy) * (1 - wz), w101 = wx * (1 - wy) * wz;
    const double w110 = wx * wy * (1 - wz), w111 = wx * wy * wz;
    return w000 * field(i, j, k, 0) + w001 * field(i, j, k + 1, 0) + w010 * field(i, j + 1, k, 0) +
           w011 * field(i, j + 1, k + 1, 0) + w100 * field(i + 1, j, k, 0) + w101 * field(i + 1, j, k + 1, 0) +
           w110 * field(i + 1, j + 1, k, 0) + w111 * field(i + 1, j + 1, k + 1, 0);
}
}   // namespace

void DiagnosticMeasurements::require_requests_open() const {
    if (ready_)
        throw std::logic_error("Diagnostic measurement requests must precede calculate()");
}
void DiagnosticMeasurements::require_ready() const {
    if (!ready_)
        throw std::logic_error("Diagnostic measurements must be calculated before writing");
}
void DiagnosticMeasurements::request_second_moment(const std::string& name, const SliceIndices& planes) {
    require_requests_open();
    include_slice_indices(moment_requests_[name], planes);
}
void DiagnosticMeasurements::request_spectrum(const std::string& name) {
    require_requests_open();
    spectrum_requests_.insert(name);
}
void DiagnosticMeasurements::request_total(TotalGridQuantity quantity, const SliceIndices& planes) {
    require_requests_open();
    include_slice_indices(total_requests_[static_cast<size_t>(quantity)], planes);
}
void DiagnosticMeasurements::request_probe(const Vector3R& point, const std::vector<std::string>& names) {
    require_requests_open();
    auto& requested = probe_requests_[{point.x(), point.y(), point.z()}];
    requested.insert(names.begin(), names.end());
}

void DiagnosticMeasurements::calculate() {
    if (ready_)
        return;
    Field3d moment;   // Reused for six components/all species; only requested planes are retained.
    for (const auto& [name, planes] : moment_requests_) {
        const auto& species = *sources_.species.at(name);
        auto& result = moments_[name];
        if (species.is_neutral() || !has_planes(planes))
            continue;
        if (moment.capacity() == 0)
            moment.resize(sources_.domain.size(), 1);
        const auto component = [&](int index, auto v1, auto v2) {
            species.calculate_pressure_component(moment, v1, v2);
            result[index] = make_field_slices(moment.sizes(), 1, planes);
            add_field_to_slices(result[index], moment);
        };
        component(0, radial_velocity, radial_velocity);
        component(1, phi_velocity, phi_velocity);
        component(2, z_velocity, z_velocity);
        component(3, radial_velocity, phi_velocity);
        component(4, radial_velocity, z_velocity);
        component(5, z_velocity, phi_velocity);
    }
    for (const auto& name : spectrum_requests_)
        spectra_.emplace(name, sources_.species.at(name)->calculate_energy_spectrum());
    for (size_t index = 0; index < total_requests_.size(); ++index) {
        if (!has_planes(total_requests_[index]))
            continue;
        const bool charge = index == static_cast<size_t>(TotalGridQuantity::ChargeDensity);
        auto& result = totals_[index];
        result = make_field_slices(sources_.domain.size(), charge ? 1 : 3, total_requests_[index]);
        for (const auto& kv : sources_.species) {
            const auto& species = *kv.second;
            if (!species.is_neutral())
                add_field_to_slices(result, charge ? species.densityOnGrid : species.currentOnGrid);
        }
    }
    for (const auto& [position, names] : probe_requests_) {
        const Vector3R point(position[0], position[1], position[2]);
        const auto normalized = sources_.domain.to_cell_coordinates(point);
        const auto en = interpolateE_linear(sources_.electric, normalized);
        const auto bn = interpolateB_linear(sources_.magnetic, normalized);
        const auto external_e = interpolateE_linear(sources_.external_electric, normalized);
        const auto external_b = interpolateB_linear(sources_.external_magnetic, normalized);
        auto& result = probes_[position];
        result.fields = {en + external_e, bn + external_b, en, bn, external_e, external_b};
        for (const auto& name : names)
            result.densities.emplace(name,
                                     sample_density(sources_.species.at(name)->densityOnGrid, point, sources_.domain));
    }
    ready_ = true;
}
const std::array<FieldSlices, 6>& DiagnosticMeasurements::second_moment(const std::string& name) const {
    require_ready();
    return moments_.at(name);
}
const EnergySpectrum& DiagnosticMeasurements::spectrum(const std::string& name) const {
    require_ready();
    return spectra_.at(name);
}
const FieldSlices& DiagnosticMeasurements::total(TotalGridQuantity quantity) const {
    require_ready();
    return totals_[static_cast<size_t>(quantity)];
}
const ProbeSample& DiagnosticMeasurements::probe(const Vector3R& point) const {
    require_ready();
    return probes_.at({point.x(), point.y(), point.z()});
}

void write_diagnostics(int timestep, const Diagnostics& values, DiagnosticSources sources,
                       const std::vector<std::unique_ptr<IDiagnosticOutput>>& outputs) {
    DiagnosticMeasurements measured(sources);
    for (const auto& output : outputs) output->request_measurements(timestep, measured);
    measured.calculate();
    for (const auto& output : outputs) output->output(timestep, values, measured);
}
