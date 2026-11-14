// ECSIM-Corr corrector passes: the velocity push is repeated with the predicted
// field moved towards (E^n+E^{n+1})/2. Each pass must shrink that mismatch and
// keep the energy balance and the charge continuity of the single-pass scheme.

#include <array>
#include <memory>
#include <stdexcept>

#include "diagnostics/Diagnostic.h"
#include "simulation/ecsim_corr/simulation_ecsim_corr.h"
#include "test_support.h"

namespace {
constexpr int STEPS = 3;

// A magnetized thermal plasma with Omega_e dt = 0.4 on cells of ten skin depths.
nlohmann::json settings(int corrector_passes) {
    nlohmann::json config = {{"Scheme", "ecsim_corr"},
                             {"Dt", 2.0},
                             {"NumCellsX", 8},
                             {"NumCellsY", 2},
                             {"NumCellsZ", 2},
                             {"Dx", 10.0},
                             {"Dy", 10.0},
                             {"Dz", 10.0},
                             {"StartFromTime", 0},
                             {"StartTimeStep", 0},
                             {"LastTimestep", STEPS},
                             {"n0", 1e10},
                             {"k_particles_reservation", 0},
                             {"Collisions", nlohmann::json::array()},
                             {"CorrectorPasses", corrector_passes}};
    config["Boundary_conditions"] = nlohmann::json::parse(
        R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])");
    config["ExternalFieldB"] = nlohmann::json::parse(R"([{"uniform_field":{"value":[0,0,0.2]}}])");
    config["diagnostics"]["outputs"] = nlohmann::json::parse(R"([{"type":"charge_continuity","interval":1}])");
    return config;
}

nlohmann::json thermal_species(const std::string& name, double charge, double mass) {
    nlohmann::json species = check::species(name, charge, mass, 1, 16);
    species["distribution"] = nlohmann::json::parse(R"([{
        "type": "initial", "sampling": "halton_quiet", "density": 1.0,
        "dist_space": {"type": "rectangle", "center": [40, 10, 10], "half_length": [40, 10, 10]},
        "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [1, 1, 1]}}])");
    return species;
}

// |(E^n+E^{n+1})/2 - E'| after the last of STEPS steps.
double field_mismatch_after_run(int corrector_passes) {
    const nlohmann::json particles = {
        {"particles", nlohmann::json::array({thermal_species("Electrons", -1, 1), thermal_species("Ions", 1, 100)})}};
    char program[] = "corrector_passes";
    char* argv[] = {program};
    std::unique_ptr<Simulation> simulation = build_simulation(settings(corrector_passes), particles, 1, argv);
    SimulationEcsimCorr& corr = dynamic_cast<SimulationEcsimCorr&>(*simulation);
    simulation->init();
    for (int step = 1; step <= STEPS; ++step) {
        simulation->prepare_step(step);
        simulation->make_step(step);
        Diagnostics values(simulation->species);
        corr.diagnostic_energy(values);
        const double kinetic = values.energy.at("Electrons") + values.energy.at("Ions");
        check::require(kinetic > 0, "the plasma must carry kinetic energy");
        check::near(values.energy.at("energyConserve") / kinetic, 0, 1e-12, "energy balance of the corrected step");
        check::require(values.charge_continuity.has_value(), "continuity was requested on every step");
        check::near(values.charge_continuity->linf, 0, 1e-12, "charge continuity of the corrected step");
    }
    const Field3d mismatch = 0.5 * (corr.fieldE + corr.fieldEn) - corr.fieldEp;
    return mismatch.norm();
}

void passes_shrink_field_mismatch() {
    std::array<double, 3> mismatch;
    for (int passes = 0; passes < 3; ++passes) mismatch[passes] = field_mismatch_after_run(passes);
    check::require(mismatch[0] > 1e-8, "the single-pass scheme must leave a predicted-field mismatch to remove");
    check::require(mismatch[1] < 0.3 * mismatch[0], "the first corrector pass must shrink the field mismatch");
    check::require(mismatch[2] < 0.3 * mismatch[1], "the second corrector pass must shrink the field mismatch");
}

void negative_pass_count_is_rejected() {
    char program[] = "corrector_passes";
    char* argv[] = {program};
    std::unique_ptr<Simulation> simulation =
        build_simulation(settings(-1), {{"particles", nlohmann::json::array()}}, 1, argv);
    bool rejected = false;
    try {
        simulation->init();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "CorrectorPasses = -1 must be rejected");
}

void scenarios() {
    check::scenario("ECSIM-Corr: corrector passes shrink the predicted-field mismatch", passes_shrink_field_mismatch);
    check::scenario("ECSIM-Corr: a negative corrector pass count is rejected", negative_pass_count_is_rejected);
}
}   // namespace

int main() {
    return check::run(scenarios);
}
