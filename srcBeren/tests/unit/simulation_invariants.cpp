#include <limits>

#include "Diagnostic.h"
#include "ParticlesDiagnostic.h"
#include "simulation_ecsim_corr.h"
#include "test_support.h"

nlohmann::json settings(const std::string& scheme) {
    return {{"Scheme", scheme},
            {"Dt", 0.25},
            {"NumCellsX", 4},
            {"NumCellsY", 4},
            {"NumCellsZ", 4},
            {"Dx", 1},
            {"Dy", 1},
            {"Dz", 1},
            {"StartFromTime", 0},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

// Each check below receives one initialized calculation. Resetting step events,
// moving particles and reading diagnostics remain explicit, in production order.
void verify_neutral_transport_and_zero_electric_charge(Simulation& simulation, ParticlesArray& gas) {
    const Vector3R start(3.95, 1.4, 1.6), velocity(0.4, 0.08, -0.04);
    gas.add_particle(Particle(start, velocity));
    gas.density_on_grid_update<ShapeType::Linear>();
    gas.kineticEnergy = gas.get_kinetic_energy();
    Field3d charge(simulation.domain.size(), 1);
    simulation.collect_charge_density(charge);
    check::near(charge.norm(), 0, 0, "neutral number density leaked into electric charge");

    simulation.prepare_step(1);
    simulation.make_step(1);
    check::require(gas.get_total_num_of_particles() == 1, "neutral lost on a periodic face");
    // dt=0.25, domain length=4: x=3.95+0.4*dt wraps to 0.05.
    for (int cell = 0; cell < gas.size(); ++cell)
        for (const auto& particle : gas.particlesData(cell)) {
            check::vector(particle.coord, Vector3R(0.05, 1.42, 1.59), 1e-13, "neutral drift and periodic wrap");
            check::vector(particle.velocity, velocity, 0, "neutral accelerated by electromagnetic fields");
        }
}

void verify_filtered_current_work_sign(SimulationEcsim& simulation, Diagnostics& diagnostic) {
    // Synthetic energy accounting: this tests the balance convention, not the
    // physical energy of a charged particle absorbed by an open wall.
    const double dt = 0.25;
    const double kinetic_before = 2, kinetic_after = 3;
    const double lost_energy = 0.5, injected_energy = 0.25;
    const double external_power = 0.2, damping_energy = 0, filtered_current_power = 0.4;
    simulation.compute_field_energy_and_conservation(
        diagnostic, simulation.bc_handler.active_range(simulation.domain.grid), dt, kinetic_before, kinetic_after,
        lost_energy, injected_energy, external_power, damping_energy, filtered_current_power);
    check::near(diagnostic.energy.at("energyResidualRaw") - diagnostic.energy.at("energyResidualDiscrete"), 0.1, 1e-14,
                "raw/discrete residual differs by dt times filtered-current work");
    check::near(diagnostic.energy.at("energyConserve"), std::abs(diagnostic.energy.at("energyResidualDiscrete")), 0,
                "legacy residual remains absolute discrete residual");
}

void verify_diagnostic_reads_preserve_step_events(SimulationEcsim& simulation, Diagnostics& diagnostic,
                                                  const ParticlesArray& gas) {
    simulation.diagnostic_energy(diagnostic);
    const auto first_read = diagnostic.energy;
    simulation.diagnostic_energy(diagnostic);
    check::require(first_read.size() == diagnostic.energy.size(), "diagnostic schema changed on repeated read");
    for (const auto& [name, value] : first_read)
        check::near(diagnostic.energy.at(name), value, 1e-15, "diagnostic read mutated step state");
    check::near(gas.diag.total_lost_energy(), 0.2, 0, "reading diagnostics consumed boundary events");
}

void neutral_transport_and_step_bookkeeping(const std::string& scheme) {
    auto config = settings(scheme);
    auto particles = nlohmann::json{{"particles", nlohmann::json::array({check::species("Gas", 0, 2)})}};
    char program[] = "simulation_invariants";
    char* argv[] = {program};
    auto simulation = build_simulation(config, particles, 1, argv);
    simulation->init();
    auto& gas = *simulation->species.at("Gas");
    auto& ecsim = dynamic_cast<SimulationEcsim&>(*simulation);

    verify_neutral_transport_and_zero_electric_charge(*simulation, gas);
    gas.diag.add_loss(Face::ZMIN, 0.2);
    Diagnostics diagnostic(simulation->species);
    ecsim.diagnostic_energy(diagnostic);
    verify_filtered_current_work_sign(ecsim, diagnostic);
    verify_diagnostic_reads_preserve_step_events(ecsim, diagnostic, gas);
    simulation->prepare_step(2);
    check::near(gas.diag.total_lost_energy(), 0, 0, "new step did not reset boundary events");
    simulation->finalize();
}
void initial_density_uses_the_scheme_charge_shape(const std::string& scheme) {
    auto cfg = check::species("Beam", -1, 1, 1, 1);
    const nlohmann::json uniform_cube = {
        {"type", "rectangle"}, {"center", {1.5, 1.5, 1.5}}, {"half_length", {0.5, 0.5, 0.5}}};
    const nlohmann::json zero_velocity = {{"type", "gaussian"}, {"mean", {0, 0, 0}}, {"sigma", {0, 0, 0}}};
    cfg["distribution"] = nlohmann::json::array(
        {{{"type", "initial"}, {"dist_space", uniform_cube}, {"dist_pulse", zero_velocity}, {"density", 1}}});
    char program[] = "initial_shape";
    char* argv[] = {program};
    auto system = settings(scheme);
    system["Boundary_conditions"] = nlohmann::json::array();
    auto sim = build_simulation(system, {{"particles", nlohmann::json::array({cfg})}}, 1, argv);
    sim->init();
    const auto& sp = *sim->species.at("Beam");
    check::require(sp.get_total_num_of_particles() == 1,
                   "initial shape regression needs a nonempty one-particle cloud");
    // Tensor-product support: Linear reaches 2^3 nodes, Quadratic 3^3.
    int support = 0;
    for (int i = 0; i < sp.densityOnGrid.capacity(); ++i) support += std::abs(sp.densityOnGrid[i]) > 0;
    check::require(support == (scheme == "ecsim_corr" ? 27 : 8),
                   "initial density uses a different shape than the scheme");
    sim->finalize();
}
void velocity_correction_adds_current_work_and_rejects_invalid_energy() {
    auto config = settings("ecsim_corr");
    char program[] = "velocity_correction";
    char* argv[] = {program};
    auto sim = build_simulation(config, {{"particles", nlohmann::json::array({check::species("Beam", -1, 1, 1, 1)})}},
                                1, argv);
    sim->init();
    auto& corr = dynamic_cast<SimulationEcsimCorr&>(*sim);
    auto& beam = *sim->species.at("Beam");
    const Vector3R velocity(1, 0, 0);
    beam.add_particle(Particle(Vector3R(2, 2, 2), velocity));
    corr.fieldE.setZero();
    corr.fieldEn.data().setConstant(2);
    corr.fieldE_external.setZero();
    beam.currentOnGrid.setZero();
    // Ebar_x = 1 and Jx = 1 in a unit volume: work over dt=1/2 is 1/2.
    const int node = GHOST_CELLS + 1;
    beam.currentOnGrid(node, node, node, Axis::X) = 1;
    const double initialEnergy = beam.get_kinetic_energy();
    corr.rescale_velocities_for_energy_balance(beam, 0.5);
    check::near(beam.get_kinetic_energy() - initialEnergy, 0.5, 1e-14,
                "velocity correction must add the deposited electric work");
    for (int cell = 0; cell < beam.size(); ++cell)
        for (Particle& particle : beam.particlesData(cell)) particle.velocity = velocity;
    // Zero/negative radicands and non-finite work must fail before any velocity changes.
    for (double current : {-2 * initialEnergy, -4 * initialEnergy, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        beam.currentOnGrid(node, node, node, Axis::X) = current;
        bool rejected = false;
        try {
            corr.rescale_velocities_for_energy_balance(beam, 0.5);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check::require(rejected, "invalid velocity correction must stop the calculation");
        for (int cell = 0; cell < beam.size(); ++cell)
            for (const Particle& particle : beam.particlesData(cell))
                check::vector(particle.velocity, velocity, 0, "failed correction changed a velocity");
    }
    sim->finalize();
}
int main() {
    return check::run([] {
        Simulation default_simulation;
        check::near(default_simulation.dt_, 0, 0, "default simulation starts with an unset timestep");
        for (const std::string scheme : {"ecsim", "ecsim_corr"}) {
            check::scenario(scheme + ": neutral transport, filtered-current work and step-event lifetime",
                            [&] { neutral_transport_and_step_bookkeeping(scheme); });
        }
        for (const std::string scheme : {"ecsim", "ecsim_corr"}) {
            check::scenario(scheme + ": initial density uses the same charge shape as the step",
                            [&] { initial_density_uses_the_scheme_charge_shape(scheme); });
        }
        check::scenario("Corr velocity correction: electric work, invalid energy and unchanged failed state",
                        velocity_correction_adds_current_work_and_rejects_invalid_energy);
    });
}
