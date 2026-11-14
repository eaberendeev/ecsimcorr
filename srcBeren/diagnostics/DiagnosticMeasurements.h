#pragma once

#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "FieldSlices.h"
#include "SpeciesMomentSource.h"
#include "particles/particle_array.h"
#include "simulation/boundary_conditions.h"

class Diagnostics;
class IDiagnosticOutput;

// Endpoint fields and particles after collision_step. Density is sampled on
// request; the deposited current retains its scheme-defined temporal level.
struct DiagnosticSources {
    const Domain& domain;
    const Species& species;
    const Field3d& electric;
    const Field3d& magnetic;
    const Field3d& external_electric;
    const Field3d& external_magnetic;
    const BoundaryConditionHandler& boundaries;
    ShapeType density_shape = SHAPE;
    // Species outside `species` (guiding centres); optional.
    const MomentSources* extra_species = nullptr;
};

enum class TotalGridQuantity { ChargeDensity, Current };
inline constexpr std::array<const char*, 6> second_moment_names{"Prr", "Ppp", "Pzz", "Prp", "Prz", "Pzp"};

// Velocity distribution of one species in a spatial box.
struct HistogramAxis {
    double min = 0.0, max = 0.0;
    int bins = 0;
};

struct DistributionRequest {
    // Kinetic: vx_vy, vz_vxy (v_z, |v_xy|), vr_vphi (cylindrical about the
    // axis of the moments); guiding centres: vpar_mu (v_par, mu = m v_perp^2 / 2|B|).
    std::string projection;
    std::array<HistogramAxis, 2> axes;
    Vector3R region_min, region_max;   // markers with region_min <= r < region_max
};

// f(a, b) = (sum of marker weights in the bin) / (da db), row-major in (a, b);
// sum f da db is the weight of the markers in the box within the axis ranges.
struct Distribution {
    DistributionRequest request;
    std::vector<double> values;
};

struct ProbeSample {
    // Total E/B, self E/B, prescribed E/B, in the existing probe column order.
    std::array<Vector3R, 6> fields;
    std::map<std::string, double> densities;
};

// One measurement frame: collect requests, calculate once, then expose const
// results to writers. Destroyed after output; no cache survives a step/restart.
class DiagnosticMeasurements {
   public:
    explicit DiagnosticMeasurements(DiagnosticSources sources) : sources_(sources) {
    }

    void request_second_moment(const std::string& species, const SliceIndices& planes);
    void request_density(const std::string& species);
    void request_spectrum(const std::string& species);
    void request_total(TotalGridQuantity quantity, const SliceIndices& planes);
    void request_probe(const Vector3R& point, const std::vector<std::string>& species);
    // `key` distinguishes several requests for one species.
    void request_distribution(const std::string& species, const std::string& key, const DistributionRequest& request);
    void calculate();

    const std::array<FieldSlices, 6>& second_moment(const std::string& species) const;
    const EnergySpectrum& spectrum(const std::string& species) const;
    const Field3d& density(const std::string& species) const;
    // Current density of a species (kinetic or extra), not a measurement.
    const Field3d& current(const std::string& species) const;
    const FieldSlices& total(TotalGridQuantity quantity) const;
    const ProbeSample& probe(const Vector3R& point) const;
    const Distribution& distribution(const std::string& species, const std::string& key) const;

   private:
    using Point = std::array<double, 3>;
    DiagnosticSources sources_;
    bool ready_ = false;
    std::map<std::string, SliceIndices> moment_requests_;
    std::set<std::string> spectrum_requests_;
    std::set<std::string> density_requests_;
    std::array<SliceIndices, 2> total_requests_;
    std::map<Point, std::set<std::string>> probe_requests_;
    std::map<std::string, std::array<FieldSlices, 6>> moments_;
    std::map<std::string, EnergySpectrum> spectra_;
    std::map<std::string, Field3d> densities_;
    std::array<FieldSlices, 2> totals_;
    std::map<Point, ProbeSample> probes_;
    std::map<std::pair<std::string, std::string>, Distribution> distributions_;   // requests, then results
    void require_requests_open() const;
    void require_ready() const;
    const SpeciesMomentSource* extra(const std::string& species) const;
    void calculate_distribution(const std::string& species, Distribution& result) const;
};

void write_diagnostics(int timestep, const Diagnostics& values, DiagnosticSources sources,
                       const std::vector<std::unique_ptr<IDiagnosticOutput>>& outputs);
