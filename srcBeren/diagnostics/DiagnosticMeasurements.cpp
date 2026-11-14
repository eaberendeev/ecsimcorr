#include "DiagnosticMeasurements.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "DensityMeasurements.h"
#include "IDiagnosticOutput.h"
#include "coupling/gather.h"
#include "coupling/moments.h"
#include "runtime/profiler.h"

namespace {
bool has_planes(const SliceIndices& planes) {
    return !planes[0].empty() || !planes[1].empty() || !planes[2].empty();
}

double sample_density(const Field3d& field, const Vector3R& point, const Domain& domain) {
    const Stencil3D<2> stencil = make_point_stencil<ShapeType::Linear>(domain.to_cell_coordinates(point));
    const Vector3I origin = to_storage_index(stencil.origin);
    const int i = origin.x(), j = origin.y(), k = origin.z();
    const double x0 = stencil(0, Axis::X), x1 = stencil(1, Axis::X);
    const double y0 = stencil(0, Axis::Y), y1 = stencil(1, Axis::Y);
    const double z0 = stencil(0, Axis::Z), z1 = stencil(1, Axis::Z);
    const double w000 = x0 * y0 * z0, w001 = x0 * y0 * z1;
    const double w010 = x0 * y1 * z0, w011 = x0 * y1 * z1;
    const double w100 = x1 * y0 * z0, w101 = x1 * y0 * z1;
    const double w110 = x1 * y1 * z0, w111 = x1 * y1 * z1;
    return w000 * field(i, j, k, 0) + w001 * field(i, j, k + 1, 0) + w010 * field(i, j + 1, k, 0) +
           w011 * field(i, j + 1, k + 1, 0) + w100 * field(i + 1, j, k, 0) + w101 * field(i + 1, j, k + 1, 0) +
           w110 * field(i + 1, j + 1, k, 0) + w111 * field(i + 1, j + 1, k + 1, 0);
}
}   // namespace

const SpeciesMomentSource* DiagnosticMeasurements::extra(const std::string& name) const {
    if (sources_.species.contains(name) || !sources_.extra_species)
        return nullptr;
    const auto it = sources_.extra_species->find(name);
    return it == sources_.extra_species->end() ? nullptr : it->second;
}

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

void DiagnosticMeasurements::request_density(const std::string& name) {
    require_requests_open();
    density_requests_.insert(name);
}

void DiagnosticMeasurements::request_total(TotalGridQuantity quantity, const SliceIndices& planes) {
    require_requests_open();
    include_slice_indices(total_requests_[static_cast<size_t>(quantity)], planes);
    if (quantity == TotalGridQuantity::ChargeDensity && has_planes(planes)) {
        for (const auto& [name, species] : sources_.species)
            if (!species->is_neutral())
                request_density(name);
        if (sources_.extra_species)
            for (const auto& [name, source] : *sources_.extra_species)
                if (!source->is_neutral())
                    request_density(name);
    }
}

void DiagnosticMeasurements::request_probe(const Vector3R& point, const std::vector<std::string>& names) {
    require_requests_open();
    auto& requested = probe_requests_[{point.x(), point.y(), point.z()}];
    requested.insert(names.begin(), names.end());
    for (const auto& name : names) request_density(name);
}

void DiagnosticMeasurements::request_distribution(const std::string& species, const std::string& key,
                                                  const DistributionRequest& request) {
    require_requests_open();
    for (const auto& axis : request.axes)
        if (!(axis.min < axis.max) || axis.bins < 1)
            throw std::invalid_argument("Distribution axes need min < max and bins >= 1");
    distributions_[{species, key}] = {request, {}};
}

namespace {
// Adds weight w at (a, b) of the histogram; outside the axis ranges is skipped.
void add_to_histogram(Distribution& result, double a, double b, double weight) {
    const auto& [x, y] = result.request.axes;
    if (!(a >= x.min && a < x.max && b >= y.min && b < y.max))
        return;
    const int i = std::min(x.bins - 1, static_cast<int>((a - x.min) / (x.max - x.min) * x.bins));
    const int j = std::min(y.bins - 1, static_cast<int>((b - y.min) / (y.max - y.min) * y.bins));
    result.values[static_cast<size_t>(i) * y.bins + j] += weight;
}

bool inside(const DistributionRequest& request, const Vector3R& r) {
    for (int d = 0; d < 3; ++d)
        if (!(r[d] >= request.region_min[d] && r[d] < request.region_max[d]))
            return false;
    return true;
}
}   // namespace

void DiagnosticMeasurements::calculate_distribution(const std::string& name, Distribution& result) const {
    const DistributionRequest& request = result.request;
    const auto& [x, y] = request.axes;
    result.values.assign(static_cast<size_t>(x.bins) * y.bins, 0.0);
    if (const SpeciesMomentSource* source = extra(name)) {
        if (request.projection != "vpar_mu")
            throw std::invalid_argument("Distribution of the guiding centres " + name + " supports vpar_mu only");
        source->for_each_marker([&](const Vector3R& r, double v_parallel, double mu, double weight) {
            if (inside(request, r))
                add_to_histogram(result, v_parallel, mu, weight);
        });
    } else {
        const ParticlesArray& species = *sources_.species.at(name);
        const Vector3R center = sources_.domain.geom.use_cylinder
                                    ? sources_.domain.geom.cyl_center
                                    : 0.5 * (sources_.domain.geom.box_min + sources_.domain.geom.box_max);
        const int projection = request.projection == "vx_vy"     ? 0
                               : request.projection == "vz_vxy"  ? 1
                               : request.projection == "vr_vphi" ? 2
                                                                 : -1;
        if (projection < 0)
            throw std::invalid_argument("Distribution of the kinetic species " + name +
                                        " supports vx_vy, vz_vxy and vr_vphi");
        const double weight = species.mpw();
        for (int cell = 0; cell < species.size(); ++cell)
            for (const Particle& p : species(cell)) {
                if (!inside(request, p.coord))
                    continue;
                const Vector3R& v = p.velocity;
                if (projection == 0) {
                    add_to_histogram(result, v.x(), v.y(), weight);
                } else if (projection == 1) {
                    add_to_histogram(result, v.z(), std::sqrt(v.x() * v.x() + v.y() * v.y()), weight);
                } else {
                    const Vector3R c = coupling::cylindrical_velocity(p.coord, v, center);
                    add_to_histogram(result, c.x(), c.y(), weight);
                }
            }
    }
    const double area = (x.max - x.min) / x.bins * (y.max - y.min) / y.bins;
    for (double& value : result.values) value /= area;
}

void DiagnosticMeasurements::calculate() {
    RECORD_TIMER;
    if (ready_)
        return;
    for (const auto& name : density_requests_) {
        if (const SpeciesMomentSource* source = extra(name)) {
            Field3d density(sources_.domain.size(), 1);
            source->deposit_density(density, sources_.density_shape);
            sources_.boundaries.apply_to_fields(density, FieldType::DENSITY, sources_.domain);
            densities_.emplace(name, std::move(density));
            continue;
        }
        densities_.emplace(name, measure_density(*sources_.species.at(name), sources_.domain, sources_.boundaries,
                                                 sources_.density_shape));
    }
    Field3d moment;   // Reused for six components/all species; only requested planes are retained.
    for (const auto& [name, planes] : moment_requests_) {
        std::array<FieldSlices, 6>& result = moments_[name];
        if (const SpeciesMomentSource* source = extra(name)) {
            if (source->is_neutral() || !has_planes(planes))
                continue;
            if (moment.capacity() == 0)
                moment.resize(sources_.domain.size(), 1);
            // Same component order as second_moment_names: rr, pp, zz, rp, rz, zp.
            constexpr int pairs[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {2, 1}};
            for (int index = 0; index < 6; ++index) {
                source->deposit_second_moment(moment, pairs[index][0], pairs[index][1]);
                result[index] = make_field_slices(moment.sizes(), 1, planes);
                add_field_to_slices(result[index], moment);
            }
            continue;
        }
        const ParticlesArray& species = *sources_.species.at(name);
        if (species.is_neutral() || !has_planes(planes))
            continue;
        if (moment.capacity() == 0)
            moment.resize(sources_.domain.size(), 1);
        const auto component = [&](int index, auto v1, auto v2) {
            coupling::deposit_second_moment(species, moment, v1, v2);
            result[index] = make_field_slices(moment.sizes(), 1, planes);
            add_field_to_slices(result[index], moment);
        };
        component(0, coupling::radial_velocity, coupling::radial_velocity);
        component(1, coupling::phi_velocity, coupling::phi_velocity);
        component(2, coupling::z_velocity, coupling::z_velocity);
        component(3, coupling::radial_velocity, coupling::phi_velocity);
        component(4, coupling::radial_velocity, coupling::z_velocity);
        component(5, coupling::z_velocity, coupling::phi_velocity);
    }
    for (const auto& name : spectrum_requests_) {
        if (extra(name))
            throw std::invalid_argument("Energy spectrum is not defined for the guiding-centre species " + name);
        spectra_.emplace(name, sources_.species.at(name)->calculate_energy_spectrum());
    }
    for (size_t index = 0; index < total_requests_.size(); ++index) {
        if (!has_planes(total_requests_[index]))
            continue;
        const bool charge = index == static_cast<size_t>(TotalGridQuantity::ChargeDensity);
        FieldSlices& result = totals_[index];
        result = make_field_slices(sources_.domain.size(), charge ? 1 : 3, total_requests_[index]);
        for (const auto& kv : sources_.species) {
            const ParticlesArray& species = *kv.second;
            if (!species.is_neutral())
                add_field_to_slices(result, charge ? densities_.at(kv.first) : species.currentOnGrid);
        }
        if (sources_.extra_species)
            for (const auto& [name, source] : *sources_.extra_species)
                if (!source->is_neutral())
                    add_field_to_slices(result, charge ? densities_.at(name) : source->current());
    }
    for (auto& [key, result] : distributions_) calculate_distribution(key.first, result);
    if (!probe_requests_.empty()) {
        require_gather_grid(sources_.electric, sources_.domain.grid);
        require_gather_grid(sources_.magnetic, sources_.domain.grid);
        require_gather_grid(sources_.external_electric, sources_.domain.grid);
        require_gather_grid(sources_.external_magnetic, sources_.domain.grid);
    }
    for (const auto& [position, names] : probe_requests_) {
        const Vector3R point(position[0], position[1], position[2]);
        const Vector3R normalized = sources_.domain.to_cell_coordinates(point);
        const Stencil3D<2> no = make_point_stencil<ShapeType::Linear>(normalized);
        const Stencil3D<2> sh = make_point_stencil<ShapeType::Linear>(normalized, Vector3R(0.5));
        const Vector3R en = interpolateE(sources_.electric, no, sh);
        const Vector3R bn = interpolateB(sources_.magnetic, no, sh);
        const Vector3R external_e = interpolateE(sources_.external_electric, no, sh);
        const Vector3R external_b = interpolateB(sources_.external_magnetic, no, sh);
        ProbeSample& result = probes_[position];
        result.fields = {en + external_e, bn + external_b, en, bn, external_e, external_b};
        for (const auto& name : names)
            result.densities.emplace(name, sample_density(densities_.at(name), point, sources_.domain));
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

const Field3d& DiagnosticMeasurements::density(const std::string& name) const {
    require_ready();
    return densities_.at(name);
}

const Field3d& DiagnosticMeasurements::current(const std::string& name) const {
    if (const SpeciesMomentSource* source = extra(name))
        return source->current();
    return sources_.species.at(name)->currentOnGrid;
}

const Distribution& DiagnosticMeasurements::distribution(const std::string& species, const std::string& key) const {
    require_ready();
    return distributions_.at({species, key});
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
    RECORD_TIMER;
    DiagnosticMeasurements measured(sources);
    for (const auto& output : outputs) output->request_measurements(timestep, measured);
    measured.calculate();
    for (const auto& output : outputs) output->output(timestep, values, measured);
}
