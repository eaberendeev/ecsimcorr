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

void exercise(const std::string& scheme) {
    auto config = settings(scheme);
    auto particles = nlohmann::json{{"particles", nlohmann::json::array({check::species("Gas", 0, 2)})}};
    char program[] = "simulation_invariants";
    char* argv[] = {program};
    auto sim = build_simulation(config, particles, 1, argv);
    sim->init();
    auto& gas = *sim->species.at("Gas");
    const Vector3R start(3.95, 1.4, 1.6), velocity(0.4, 0.08, -0.04);
    gas.add_particle(Particle(start, velocity));
    gas.density_on_grid_update<ShapeType::Linear>();
    gas.kineticEnergy = gas.get_kinetic_energy();
    Field3d rho(sim->domain.size(), 1);
    sim->collect_charge_density(rho);
    check::near(rho.norm(), 0, 0, "neutral number density leaked into electric charge");
    sim->prepare_step(1);
    sim->make_step(1);
    check::require(gas.get_total_num_of_particles() == 1, "neutral lost on a periodic face");
    for (int cell = 0; cell < gas.size(); ++cell)
        for (const auto& p : gas.particlesData(cell)) {
            check::vector(p.coord, Vector3R(0.05, 1.42, 1.59), 1e-13, "neutral drift and periodic wrap");
            check::vector(p.velocity, velocity, 0, "neutral accelerated by electromagnetic fields");
        }
    auto& ecsim = dynamic_cast<SimulationEcsim&>(*sim);
    gas.diag.add_loss(Face::ZMIN, 0.2);
    Diagnostics diag(sim->species);
    ecsim.diagnostic_energy(diag);
    // A synthetic filtered-current work checks the sign and preserves the old
    // absolute diagnostic independently of the charged-particle dynamics.
    ecsim.compute_field_energy_and_conservation(diag, ecsim.bc_handler.active_range(sim->domain.grid), 0.25, 2, 3, 0.5,
                                                0.25, 0.2, 0, 0.4);
    check::near(diag.energy.at("energyResidualRaw") - diag.energy.at("energyResidualDiscrete"), 0.1, 1e-14,
                "raw/discrete residual differs by dt times filtered-current work");
    check::near(diag.energy.at("energyConserve"), std::abs(diag.energy.at("energyResidualDiscrete")), 0,
                "legacy residual remains absolute discrete residual");
    ecsim.diagnostic_energy(diag);
    const auto energy = diag.energy;
    ecsim.diagnostic_energy(diag);
    check::require(energy.size() == diag.energy.size(), "diagnostic schema changed on repeated read");
    for (const auto& item : energy)
        check::near(diag.energy.at(item.first), item.second, 1e-15, "diagnostic read mutated step state");
    check::near(gas.diag.total_lost_energy(), 0.2, 0, "reading diagnostics consumed boundary events");
    sim->prepare_step(2);
    check::near(gas.diag.total_lost_energy(), 0, 0, "new step did not reset boundary events");
    sim->finalize();
}
void initial_shape(const std::string& scheme) {
    auto cfg = check::species("Beam", -1, 1, 1, 1);
    cfg["distribution"] = nlohmann::json::parse(
        R"([{"type":"initial","dist_space":{"type":"rectangle","center":[1.5,1.5,1.5],"half_length":[0.5,0.5,0.5]},"dist_pulse":{"type":"gaussian","mean":[0,0,0],"sigma":[0,0,0]},"density":1}])");
    char program[] = "initial_shape";
    char* argv[] = {program};
    auto system = settings(scheme);
    system["Boundary_conditions"] = nlohmann::json::array();
    auto sim = build_simulation(system, {{"particles", nlohmann::json::array({cfg})}}, 1, argv);
    sim->init();
    const auto& sp = *sim->species.at("Beam");
    check::require(sp.get_total_num_of_particles() == 1,
                   "initial shape regression needs a nonempty one-particle cloud");
    int support = 0;
    for (int i = 0; i < sp.densityOnGrid.capacity(); ++i) support += std::abs(sp.densityOnGrid[i]) > 0;
    check::require(support == (scheme == "ecsim_corr" ? 27 : 8),
                   "initial density uses a different shape than the scheme");
    sim->finalize();
}
void velocity_correction() {
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
    corr.correctv(beam, 0.5);
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
            corr.correctv(beam, 0.5);
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
        exercise("ecsim");
        exercise("ecsim_corr");
        initial_shape("ecsim");
        initial_shape("ecsim_corr");
        velocity_correction();
    });
}
