#pragma once

#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "FieldSlices.h"
#include "ParticlesArray.h"

class Diagnostics;
class IDiagnosticOutput;

// Endpoint fields and particles after collision_step. Deposited rho/J retain
// their scheme-defined levels; measuring does not redeposit or move particles.
struct DiagnosticSources {
    const Domain& domain;
    const Species& species;
    const Field3d& electric;
    const Field3d& magnetic;
    const Field3d& external_electric;
    const Field3d& external_magnetic;
};

enum class TotalGridQuantity { ChargeDensity, Current };
inline constexpr std::array<const char*, 6> second_moment_names{"Prr", "Ppp", "Pzz", "Prp", "Prz", "Pzp"};
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
    void request_spectrum(const std::string& species);
    void request_total(TotalGridQuantity quantity, const SliceIndices& planes);
    void request_probe(const Vector3R& point, const std::vector<std::string>& species);
    void calculate();

    const std::array<FieldSlices, 6>& second_moment(const std::string& species) const;
    const EnergySpectrum& spectrum(const std::string& species) const;
    const FieldSlices& total(TotalGridQuantity quantity) const;
    const ProbeSample& probe(const Vector3R& point) const;

   private:
    using Point = std::array<double, 3>;
    DiagnosticSources sources_;
    bool ready_ = false;
    std::map<std::string, SliceIndices> moment_requests_;
    std::set<std::string> spectrum_requests_;
    std::array<SliceIndices, 2> total_requests_;
    std::map<Point, std::set<std::string>> probe_requests_;
    std::map<std::string, std::array<FieldSlices, 6>> moments_;
    std::map<std::string, EnergySpectrum> spectra_;
    std::array<FieldSlices, 2> totals_;
    std::map<Point, ProbeSample> probes_;
    void require_requests_open() const;
    void require_ready() const;
};

void write_diagnostics(int timestep, const Diagnostics& values, DiagnosticSources sources,
                       const std::vector<std::unique_ptr<IDiagnosticOutput>>& outputs);
