// ECSIM-Path: ECSIM-Corr without the energy correction. Each particle is pushed
// with the mean field on its straight step path and deposits that path; the
// missing magnetic moment of the path is a magnetization current. The step
// energy balance is not exact (the predictor freezes the path weights), but it
// must stay small, the magnetization exchange must be exact, charge continuity
// exact, and the magnetization must vanish without a magnetic field.

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "diagnostics/Diagnostic.h"
#include "simulation/ecsim_path/simulation_ecsim_path.h"
#include "test_support.h"

namespace {
constexpr int STEPS = 3;

nlohmann::json settings(double field, double dt) {
    nlohmann::json config = {{"Scheme", "ecsim_path"},
                             {"Dt", dt},
                             {"NumCellsX", 8},
                             {"NumCellsY", 4},
                             {"NumCellsZ", 4},
                             {"Dx", 1.0},
                             {"Dy", 1.0},
                             {"Dz", 1.0},
                             {"StartFromTime", 0},
                             {"StartTimeStep", 0},
                             {"LastTimestep", STEPS},
                             {"n0", 1e10},
                             {"k_particles_reservation", 0},
                             {"Collisions", nlohmann::json::array()}};
    config["Boundary_conditions"] = nlohmann::json::parse(
        R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])");
    config["ExternalFieldB"] = {{{"uniform_field", {{"value", {0.0, 0.3 * field, field}}}}}};
    config["diagnostics"]["outputs"] = nlohmann::json::parse(R"([{"type":"charge_continuity","interval":1}])");
    return config;
}

// 5 keV: a thermal path of 0.1 cell per unit time.
nlohmann::json thermal_species(const std::string& name, double charge, double mass) {
    nlohmann::json species = check::species(name, charge, mass, 1, 16);
    species["distribution"] = nlohmann::json::parse(R"([{
        "type": "initial", "sampling": "halton_quiet", "density": 1.0,
        "dist_space": {"type": "rectangle", "center": [4, 2, 2], "half_length": [4, 2, 2]},
        "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [5, 5, 5]}}])");
    return species;
}

struct StepValues {
    double kinetic = 0, balance = 0, magnetization = 0, remainder = 0, continuity = 0;
};

StepValues worst_step(double field, double dt) {
    const nlohmann::json particles = {
        {"particles", nlohmann::json::array({thermal_species("Electrons", -1, 1), thermal_species("Ions", 1, 100)})}};
    char program[] = "ecsim_path";
    char* argv[] = {program};
    std::unique_ptr<Simulation> simulation = build_simulation(settings(field, dt), particles, 1, argv);
    SimulationEcsimPath& path = dynamic_cast<SimulationEcsimPath&>(*simulation);
    simulation->init();
    StepValues worst;
    for (int step = 1; step <= STEPS; ++step) {
        simulation->prepare_step(step);
        simulation->make_step(step);
        Diagnostics values(simulation->species);
        path.diagnostic_values(values);
        const double kinetic = values.values.at("Electrons") + values.values.at("Ions");
        check::require(kinetic > 0, "the plasma must carry kinetic energy");
        check::require(values.charge_continuity.has_value(), "continuity was requested on every step");
        worst.kinetic = kinetic;
        worst.balance = std::max(worst.balance, std::abs(values.values.at("energyConserve")) / kinetic);
        worst.magnetization = std::max(worst.magnetization, std::abs(values.values.at("magnetizationWork")) / kinetic);
        worst.remainder = std::max(worst.remainder, std::abs(values.values.at("magnetizationRemainder")) / kinetic);
        worst.continuity = std::max(worst.continuity, values.charge_continuity->linf);
    }
    return worst;
}

void unmagnetized_plasma() {
    const StepValues worst = worst_step(0.0, 2.0);
    check::near(worst.continuity, 0, 1e-12, "charge continuity of the path current");
    check::equal(worst.magnetization, 0.0, 0.0, "no magnetization without a magnetic field");
    check::require(worst.balance < 1e-6, "small step energy error, got " + std::to_string(worst.balance));
}

// Omega_e dt = 10: the straight path carries 4 % of the magnetic moment; the
// rest is the magnetization current, whose exchange must keep the balance.
void strongly_magnetized_plasma() {
    const StepValues worst = worst_step(5.0, 2.0);
    check::near(worst.continuity, 0, 1e-12, "charge continuity of the path current");
    check::require(worst.magnetization > 0, "the magnetization exchanges energy");
    check::require(worst.balance < 1e-6, "small step energy error, got " + std::to_string(worst.balance));
    check::require(worst.balance < 1e-2 * worst.magnetization,
                   "the magnetization exchange is exact: the balance error is far below its work");
}

void invalid_settings_are_rejected() {
    char program[] = "ecsim_path";
    char* argv[] = {program};
    for (const std::string key : {"CorrectorPasses", "EnergyCorrection"}) {
        nlohmann::json config = settings(0.2, 2.0);
        config[key] = key == "CorrectorPasses" ? nlohmann::json(0) : nlohmann::json("global");
        std::unique_ptr<Simulation> simulation =
            build_simulation(config, {{"particles", nlohmann::json::array()}}, 1, argv);
        bool rejected = false;
        try {
            simulation->init();
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check::require(rejected, key + " must be rejected");
    }
}

void scenarios() {
    check::scenario("ECSIM-Path: unmagnetized plasma, no magnetization, small energy error",
                    unmagnetized_plasma);
    check::scenario("ECSIM-Path: Omega_e dt = 10, the magnetization exchange keeps the balance",
                    strongly_magnetized_plasma);
    check::scenario("ECSIM-Path: no pass and an energy correction are rejected", invalid_settings_are_rejected);
}
}   // namespace

int main() {
    return check::run(scenarios);
}
