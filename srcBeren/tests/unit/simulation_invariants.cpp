#include <limits>
#include <type_traits>

#include "diagnostics/DensityMeasurements.h"
#include "diagnostics/Diagnostic.h"
#include "diagnostics/EnergyMeasurements.h"
#include "simulation/ecsim/simulation_ecsim.h"
#include "simulation/ecsim_corr/simulation_ecsim_corr.h"
#include "test_support.h"

static_assert(std::is_abstract_v<SimulationEcsimBase>);
static_assert(!std::is_base_of_v<SimulationEcsim, SimulationEcsimCorr>);
static_assert(!std::is_base_of_v<SimulationEcsimCorr, SimulationEcsim>);

// A borrowed domain must never bind to a temporary owner.
static_assert(!std::is_constructible_v<ParticlesArray, const nlohmann::json&, Domain&&>);
static_assert(!std::is_constructible_v<ParticlesArray, const nlohmann::json&, const Domain&&>);
// Species borrow Simulation::domain; the owner must keep its address throughout their lifetime.
static_assert(!std::is_copy_constructible_v<Simulation>);
static_assert(!std::is_move_constructible_v<Simulation>);

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

void verify_one_particle_in_its_cell(const ParticlesArray& species, const Vector3R& position, const Vector3R& velocity,
                                     const std::string& message) {
    check::equal(species.get_total_num_of_particles(), 1, message + ": particle count");
    const auto& domain = species.get_domain();
    const auto expected_cell = domain.storage_cell_index(position);
    const int expected_index = domain.sind(expected_cell.x(), expected_cell.y(), expected_cell.z());
    for (int cell = 0; cell < species.size(); ++cell)
        for (const auto& particle : species(cell)) {
            check::equal(cell, expected_index, message + ": particle storage follows its coordinate");
            check::vector(particle.coord, position, 0, message + ": position");
            check::vector(particle.velocity, velocity, 0, message + ": velocity");
        }
}

void full_steps_drift_neutrals_once_with_disabled_reactions(const std::string& scheme, int channel_count) {
    auto config = settings(scheme);
    for (int channel = 0; channel < channel_count; ++channel)
        config["Collisions"].push_back({{"type", "neutral_ionization"},
                                        {"charged_species", channel == 0 ? "Projectile" : "Countercharge"},
                                        {"neutral_species", "Gas"},
                                        {"electron_product", "SecondaryElectrons"},
                                        {"ion_product", "ProductIons"},
                                        {"scheme", "physical_only"},
                                        {"electron_ionization", false},
                                        {"proton_ionization", false},
                                        {"proton_charge_exchange", false}});
    const auto particles =
        nlohmann::json{{"particles", nlohmann::json::array({check::species("Gas", 0, 1837, 1, 1),
                                                            check::species("Projectile", -1, 1, 1, 1),
                                                            check::species("Countercharge", 1, 1836, 1, 1),
                                                            check::species("SecondaryElectrons", -1, 1, 1, 1),
                                                            check::species("ProductIons", 1, 1836, 1, 1)})}};
    char program[] = "neutral_full_steps";
    char* argv[] = {program};
    auto simulation = build_simulation(config, particles, 1, argv);
    simulation->init();
    check::equal(simulation->collision_manager.enabled(), channel_count > 0,
                 "configured neutral channels enable the collision stage");
    auto& gas = *simulation->species.at("Gas");
    auto& projectile = *simulation->species.at("Projectile");
    auto& countercharge = *simulation->species.at("Countercharge");
    auto& ecsim = dynamic_cast<SimulationEcsimBase&>(*simulation);
    // The stationary opposite charges cancel rho and J exactly. After the
    // first periodic crossing Gas shares their cell, exercising both disabled
    // charged-neutral reaction loops without electromagnetic acceleration.
    const Vector3R charged_position(0.5, 1.5, 1.5), stationary(0, 0, 0);
    projectile.add_particle(Particle(charged_position, stationary));
    countercharge.add_particle(Particle(charged_position, stationary));
    const Vector3R start(3.9375, 0.9375, 1.5), velocity(0.5, 0.25, -0.125);
    gas.add_particle(Particle(start, velocity));
    // Dyadic positions, dt and velocities make r_n=r_0+n*dt*v exact.
    // One step crosses both the physical XMAX face and an interior Y cell face;
    // each displacement is smaller than one cell, as update_cells requires.
    for (int step = 1; step <= 4; ++step) {
        Vector3R expected = start + step * simulation->dt_ * velocity;
        for (int axis = 0; axis < 3; ++axis) {
            const double lower = simulation->domain.geom.box_min[axis];
            const double length = simulation->domain.geom.box_max[axis] - lower;
            expected[axis] -= length * std::floor((expected[axis] - lower) / length);
        }
        const std::string label = "step " + std::to_string(step);
        simulation->prepare_step(step);
        simulation->make_step(step);
        verify_one_particle_in_its_cell(gas, expected, velocity, label + " after the scheme step");
        const auto gas_cell = simulation->domain.storage_cell_index(expected);
        const int cell = simulation->domain.sind(gas_cell.x(), gas_cell.y(), gas_cell.z());
        check::equal(projectile(cell).size(), size_t{1}, label + ": gas shares the electron projectile cell");
        check::equal(countercharge(cell).size(), size_t{1}, label + ": gas shares the ion projectile cell");

        simulation->collision_step(step);
        verify_one_particle_in_its_cell(gas, expected, velocity, label + " after the collision stage");
        verify_one_particle_in_its_cell(projectile, charged_position, stationary, label + " electron projectile");
        verify_one_particle_in_its_cell(countercharge, charged_position, stationary, label + " ion projectile");
        for (const std::string product : {"SecondaryElectrons", "ProductIons"})
            check::equal(simulation->species.at(product)->get_total_num_of_particles(), 0,
                         label + ": disabled reactions produced " + product);
        const auto charge = measure_charge_density(simulation->species, simulation->domain, simulation->bc_handler,
                                                   simulation->diagnostic_density_shape());
        check::near(charge.norm(), 0, 0, label + ": neutral number density does not contribute to electric charge");
        check::near(gas.currentOnGrid.norm(), 0, 0, label + ": neutral drift produces no electric current");
        check::near(ecsim.fieldJp.norm(), 0, 0, label + ": stationary charges and gas produce no predictor current");
        if (scheme == "ecsim_corr") {
            const auto& corr = dynamic_cast<const SimulationEcsimCorr&>(*simulation);
            check::near(corr.fieldJe.norm(), 0, 0,
                        label + ": stationary charges and gas produce no correction current");
        }
    }
    simulation->finalize();
}

// Each check below receives one initialized calculation. Resetting step events,
// moving particles and reading diagnostics remain explicit, in production order.
void verify_neutral_transport_and_zero_electric_charge(Simulation& simulation, ParticlesArray& gas) {
    const Vector3R start(3.95, 1.4, 1.6), velocity(0.4, 0.08, -0.04);
    gas.add_particle(Particle(start, velocity));
    gas.step_state.kinetic_energy_before_step = gas.get_kinetic_energy();
    const auto charge = measure_charge_density(simulation.species, simulation.domain, simulation.bc_handler,
                                               simulation.diagnostic_density_shape());
    check::near(charge.norm(), 0, 0, "neutral number density leaked into electric charge");

    simulation.prepare_step(1);
    simulation.make_step(1);
    check::require(gas.get_total_num_of_particles() == 1, "neutral lost on a periodic face");
    // dt=0.25, domain length=4: x=3.95+0.4*dt wraps to 0.05.
    for (int cell = 0; cell < gas.size(); ++cell)
        for (const Particle& particle : gas(cell)) {
            check::vector(particle.coord, Vector3R(0.05, 1.42, 1.59), 1e-13, "neutral drift and periodic wrap");
            check::vector(particle.velocity, velocity, 0, "neutral accelerated by electromagnetic fields");
        }
}

void verify_filtered_current_work_sign(SimulationEcsimBase& simulation, Diagnostics& diagnostic) {
    // Synthetic energy accounting: this tests the balance convention, not the
    // physical energy of a charged particle absorbed by an open wall.
    const double dt = 0.25;
    const double kinetic_before = 2, kinetic_after = 3;
    const double lost_energy = 0.5, injected_energy = 0.25;
    const double external_power = 0.2, damping_energy = 0, filtered_current_power = 0.4;
    const EnergyFields fields{simulation.fieldE,          simulation.fieldEn, simulation.fieldEp,
                              simulation.fieldB,          simulation.fieldBn, simulation.fieldBInit,
                              simulation.fieldE_external, simulation.fieldJp};
    measure_field_energy_and_conservation(
        diagnostic, fields, simulation.bc_handler.active_range(simulation.domain.grid), dt, kinetic_before,
        kinetic_after, lost_energy, injected_energy, external_power, damping_energy, filtered_current_power);
    check::near(diagnostic.energy.at("energyResidualRaw") - diagnostic.energy.at("energyResidualDiscrete"), 0.1, 1e-14,
                "raw/discrete residual differs by dt times filtered-current work");
    check::near(diagnostic.energy.at("energyConserve"), std::abs(diagnostic.energy.at("energyResidualDiscrete")), 0,
                "legacy residual remains absolute discrete residual");
}

void verify_diagnostic_reads_preserve_step_events(SimulationEcsimBase& simulation, Diagnostics& diagnostic,
                                                  const ParticlesArray& gas) {
    simulation.diagnostic_energy(diagnostic);
    const auto first_read = diagnostic.energy;
    simulation.diagnostic_energy(diagnostic);
    check::require(first_read.size() == diagnostic.energy.size(), "diagnostic schema changed on repeated read");
    for (const auto& [name, value] : first_read)
        check::near(diagnostic.energy.at(name), value, 1e-15, "diagnostic read mutated step state");
    check::near(gas.step_state.total_lost_energy(), 0.2, 0, "reading diagnostics consumed boundary events");
}

void neutral_transport_and_step_bookkeeping(const std::string& scheme) {
    auto config = settings(scheme);
    auto particles = nlohmann::json{{"particles", nlohmann::json::array({check::species("Gas", 0, 2)})}};
    char program[] = "simulation_invariants";
    char* argv[] = {program};
    auto simulation = build_simulation(config, particles, 1, argv);
    simulation->init();
    auto& gas = *simulation->species.at("Gas");
    check::require(&gas.get_domain() == &simulation->domain, "neutral species copied the simulation domain");
    auto& ecsim = dynamic_cast<SimulationEcsimBase&>(*simulation);

    verify_neutral_transport_and_zero_electric_charge(*simulation, gas);
    gas.step_state.add_loss(Face::ZMIN, 0.2);
    Diagnostics diagnostic(simulation->species);
    ecsim.diagnostic_energy(diagnostic);
    verify_filtered_current_work_sign(ecsim, diagnostic);
    verify_diagnostic_reads_preserve_step_events(ecsim, diagnostic, gas);
    simulation->prepare_step(2);
    check::near(gas.step_state.total_lost_energy(), 0, 0, "new step did not reset boundary events");
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
    check::require(&sp.get_domain() == &sim->domain, "charged species copied the simulation domain");
    check::require(sp.get_total_num_of_particles() == 1,
                   "initial shape regression needs a nonempty one-particle cloud");
    // Tensor-product support: Linear reaches 2^3 nodes, Quadratic 3^3.
    int support = 0;
    const auto density = measure_density(sp, sim->domain, sim->bc_handler, sim->diagnostic_density_shape());
    for (int i = 0; i < density.capacity(); ++i) support += std::abs(density[i]) > 0;
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
    beam.currentOnGrid(node, node, node, static_cast<int>(Axis::X)) = 1;
    const double initialEnergy = beam.get_kinetic_energy();
    for (double predictor_power : {0.0, 0.25, -0.25}) {
        for (int cell = 0; cell < beam.size(); ++cell)
            for (Particle& particle : beam(cell)) particle.velocity = velocity;
        corr.rescale_velocities_for_energy_balance(beam, 0.5, predictor_power);
        check::near(beam.get_kinetic_energy() - initialEnergy, 0.5 - 0.5 * predictor_power, 1e-14,
                    "velocity correction must add deposited work minus the supplied Boris work");
    }
    for (int cell = 0; cell < beam.size(); ++cell)
        for (Particle& particle : beam(cell)) particle.velocity = velocity;
    // Zero/negative radicands and non-finite work must fail before any velocity changes.
    for (double current : {-2 * initialEnergy, -4 * initialEnergy, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        beam.currentOnGrid(node, node, node, static_cast<int>(Axis::X)) = current;
        bool rejected = false;
        try {
            corr.rescale_velocities_for_energy_balance(beam, 0.5, 0);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check::require(rejected, "invalid velocity correction must stop the calculation");
        for (int cell = 0; cell < beam.size(); ++cell)
            for (const Particle& particle : beam(cell))
                check::vector(particle.velocity, velocity, 0, "failed correction changed a velocity");
    }
    sim->finalize();
}

void energy_measurements_preserve_charged_state(const std::string& scheme) {
    auto system = settings(scheme);
    char program[] = "energy_measurements";
    char* argv[] = {program};
    auto sim =
        build_simulation(system, {{"particles", nlohmann::json::array({check::species("Beam", -1, 2)})}}, 1, argv);
    sim->init();
    auto& ecsim = dynamic_cast<SimulationEcsimBase&>(*sim);
    auto& beam = *sim->species.at("Beam");
    beam.add_particle(Particle(Vector3R(2, 2, 2), Vector3R(1, 2, 3)));
    beam.currentOnGrid.data().setConstant(0.75);
    beam.step_state.add_loss(Face::ZMIN, 0.2);
    ecsim.fieldE_external.data().setConstant(0.1);
    ecsim.fieldBn.data().setConstant(0.2);
    ecsim.fieldBInit.data().setConstant(0.3);
    ecsim.fieldBFull.data().setConstant(9);
    const Field3d current = beam.currentOnGrid;
    const Field3d predictor_current = ecsim.fieldJp;
    const Field3d magnetic_scratch = ecsim.fieldBFull;
    Diagnostics diagnostic(sim->species);
    ecsim.diagnostic_energy(diagnostic);
    const auto first = diagnostic.energy;
    ecsim.diagnostic_energy(diagnostic);
    for (const auto& [name, value] : first)
        check::near(diagnostic.energy.at(name), value, 0, "repeated measurement changed " + name);
    check::near((beam.currentOnGrid - current).norm(), 0, 0, "energy measurement changed deposited current");
    check::near((ecsim.fieldJp - predictor_current).norm(), 0, 0, "energy measurement changed predictor current");
    check::near((ecsim.fieldBFull - magnetic_scratch).norm(), 0, 0, "energy measurement changed numerical scratch");
    check::near(diagnostic.energy.at("Beam"), beam.mass() * beam.mpw() * 7, 1e-13,
                "energy measurement must use the current particle velocity");
    check::near(beam.step_state.total_lost_energy(), 0.2, 0, "energy measurement consumed loss events");
    sim->finalize();
}

void damping_parameters_belong_to_the_initialized_fields(const std::string& scheme) {
    nlohmann::json config = settings(scheme);
    config["DampingType"] = "Rectangle";
    config["DampingCoeff"] = 0.8;
    config["DampCellsX_glob"] = {1, 0};
    config["DampCellsY_glob"] = {0, 0};
    config["DampCellsZ_glob"] = {0, 0};
    char program[] = "damping_lifetime";
    char* argv[] = {program};
    auto simulation = build_simulation(config, {{"particles", nlohmann::json::array()}}, 1, argv);
    simulation->init();
    SimulationEcsimBase& ecsim = dynamic_cast<SimulationEcsimBase&>(*simulation);
    ecsim.fieldEn.data().setConstant(2);
    ecsim.fieldBn.data().setConstant(3);
    simulation->system_config["DampingType"] = "unknown";
    simulation->system_config["DampingCoeff"] = 0;
    for (int step : {1, 2}) {
        simulation->prepare_step(step);
        const double factor = std::pow(0.8, step);
        check::near(ecsim.fieldE(0, 3, 3, 0), 2 * factor, 1e-14, "initialized damping applies to E endpoint");
        check::near(ecsim.fieldB(0, 3, 3, 0), 3 * factor, 1e-14, "initialized damping applies to B endpoint");
        check::near(ecsim.fieldE(3, 3, 3, 0), 2, 0, "interior E remains unchanged");
        const double face_nodes = simulation->domain.size().y() * simulation->domain.size().z();
        const double expected_loss = 0.5 * (4 + 9) * 3 * face_nodes * (1 - 0.8 * 0.8) * std::pow(0.8, 2 * (step - 1));
        check::near(ecsim.dampingEnergy_, expected_loss, 1e-10, "damping energy belongs to the current step");
    }
    simulation->finalize();
}

void injection_loss_emission_and_diagnostic_lifetime(const std::string& scheme, bool neutral, bool emission) {
    nlohmann::json system = settings(scheme);
    system["Boundary_conditions"] = nlohmann::json::parse(
        R"([{"open":[{"face":"XMIN"},{"face":"XMAX"}]},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])");
    const std::string source_name = neutral ? "Gas" : "Ions";
    if (emission)
        system["Boundary_conditions"].push_back(
            {{"second_emission",
              {{"face", "XMAX"},
               {"product", "Electrons"},
               {"sources", nlohmann::json::array({{{"species", source_name},
                                                   {"yield", 1},
                                                   {"energy", {{"type", "fraction"}, {"fraction", 0.25}}}}})}}}});
    nlohmann::json config = check::species(source_name, neutral ? 0 : 1, 4, 1e-8, 1);
    // Volume .002 and density 1000 produce exactly two macroparticles per step.
    config["distribution"] = nlohmann::json::array(
        {{{"type", "injection"},
          {"density", 1000},
          {"dist_space", {{"type", "rectangle"}, {"center", {3.99, 2, 2}}, {"half_length", {0.001, 0.5, 0.5}}}},
          {"dist_pulse", {{"type", "gaussian"}, {"mean", {0.5, 0, 0}}, {"sigma", {0, 0, 0}}}}}});
    char program[] = "source_events";
    char* argv[] = {program};
    std::unique_ptr<Simulation> simulation = build_simulation(
        system, {{"particles", nlohmann::json::array({config, check::species("Electrons", -1, 1, 1e-8, 1)})}}, 1, argv);
    simulation->init();
    ParticlesArray& source = *simulation->species.at(source_name);
    ParticlesArray& product = *simulation->species.at("Electrons");
    simulation->prepare_step(1);
    check::near(source.step_state.kinetic_energy_before_step, 0, 0, "beginning energy is captured before injection");
    check::equal(source.get_total_num_of_particles(), 2, "two injected particles");
    check::equal(source.step_state.injection_count, 2, "two injection accounting events");
    check::near(source.step_state.injection_energy, 1e-8, 1e-22, "injected energy = N*m*mpw*v^2/2");
    // Execute drift and physical wall interaction explicitly, without a field
    // solve. Exact kinematics isolate source/loss/emission bookkeeping.
    source.move(simulation->dt_);
    simulation->bc_handler.apply_to_particles(source, simulation->species, simulation->domain);
    simulation->bc_handler.flush_species(simulation->species);
    check::equal(source.get_total_num_of_particles(), 0, "all injected particles are absorbed");
    check::near(source.step_state.total_lost_energy(), 1e-8, 1e-22, "absorbed injected energy");
    check::equal(source.step_state.boundary.at(Face::XMAX).lost_count, 2, "two absorption events");
    check::equal(product.get_total_num_of_particles(), emission ? 2 : 0, "emission yield times absorbed count");
    if (emission) {
        check::near(product.step_state.boundary.at(Face::XMAX).emitted_energy, 2.5e-9, 1e-22,
                    "quarter incident energy emitted");
        check::near(product.get_kinetic_energy(), 2.5e-9, 1e-22, "emission events match product kinetic energy");
    }
    SimulationEcsimBase& ecsim = dynamic_cast<SimulationEcsimBase&>(*simulation);
    Diagnostics diagnostic(simulation->species);
    ecsim.diagnostic_energy(diagnostic);
    const std::map<std::string, double> first_read = diagnostic.energy;
    ecsim.diagnostic_energy(diagnostic);
    for (const auto& [key, value] : first_read)
        check::near(diagnostic.energy.at(key), value, 0, "reading diagnostics preserves source/loss/emission state");
    check::near(diagnostic.energy.at("energyResidualRaw"), 0, 1e-21, "injection/loss/emission energy balance");
    simulation->prepare_step(2);
    check::near(source.step_state.total_lost_energy(), 0, 0, "new step resets old losses before injection");
    check::near(product.step_state.kinetic_energy_before_step, emission ? 2.5e-9 : 0, 1e-22,
                "next step captures existing emitted particles");
    check::require(product.step_state.boundary.empty(), "new step resets old emission events");
    simulation->finalize();
void small_periodic_energy() {
    for (int cells : {1, 2}) {
        auto cfg = settings("ecsim_corr");
        cfg["NumCellsZ"] = cells;
        char program[] = "small_periodic_energy";
        char* argv[] = {program};
        auto sim = build_simulation(cfg,
                                    {{"particles", nlohmann::json::array({check::species("Electrons", -1, 1, 1, 1),
                                                                          check::species("Ions", 1, 9, 1, 1)})}},
                                    1, argv);
        sim->init();
        auto& fields = dynamic_cast<SimulationEcsimCorr&>(*sim);
        const Vector3R position(3.95, 3.96, cells - 0.03);
        sim->species.at("Electrons")->add_particle(Particle(position, {0.4, 0.08, 0.2}));
        sim->species.at("Ions")->add_particle(Particle(position, {-0.04, -0.03, -0.02}));
        const auto range = sim->bc_handler.active_range(sim->domain.grid);
        auto energy = [&] {
            double total = 0.5 * dot_product_sum(fields.fieldEn, fields.fieldEn, range) +
                           0.5 * dot_product_sum(fields.fieldBn, fields.fieldBn, range);
            for (const auto& [name, sp] : sim->species) total += sp->get_kinetic_energy();
            return total;
        };
        const double initial = energy();
        for (int step = 1; step <= 10; ++step) {
            sim->prepare_step(step);
            sim->make_step(step);
            check::near(energy(), initial, 1e-9, "ECSIM-Corr energy with one/two periodic Z cells");
            const auto charge = measure_charge_density(sim->species, sim->domain, sim->bc_handler,
                                                       sim->diagnostic_density_shape());
            Field3d gauss(sim->domain.size(), 1);
            gauss.data() = sim->maxwell.divE * fields.fieldEn.data() - charge.data();
            check::near(gauss.norm(), 0, 3e-8, "Gauss law across thin periodic boundaries");
            for (const auto& [name, sp] : sim->species)
                check::require(sp->get_total_num_of_particles() == 1, "particle lost across thin periodic seam");
        }
        sim->finalize();
    }
}

int main() {
    return check::run([] {
        Simulation default_simulation;
        check::near(default_simulation.dt_, 0, 0, "default simulation starts with an unset timestep");
        for (const std::string scheme : {"ecsim", "ecsim_corr"}) {
            for (int channels : {0, 1, 2})
                check::scenario(scheme + ": one neutral drift per full step with " + std::to_string(channels) +
                                    " disabled reaction channels",
                                [&] { full_steps_drift_neutrals_once_with_disabled_reactions(scheme, channels); });
            check::scenario(scheme + ": damping parameters are fixed at field initialization",
                            [&] { damping_parameters_belong_to_the_initialized_fields(scheme); });
            for (bool neutral : {false, true})
                check::scenario(scheme + (neutral ? ": neutral injection and loss" : ": charged injection and loss"),
                                [&] { injection_loss_emission_and_diagnostic_lifetime(scheme, neutral, false); });
            check::scenario(scheme + ": injection, absorption, secondary emission and repeated diagnostic reads",
                            [&] { injection_loss_emission_and_diagnostic_lifetime(scheme, false, true); });
            check::scenario(scheme + ": neutral transport, filtered-current work and step-event lifetime",
                            [&] { neutral_transport_and_step_bookkeeping(scheme); });
        }
        for (const std::string scheme : {"ecsim", "ecsim_corr"}) {
            check::scenario(scheme + ": read-only energy measurement of charged state",
                            [&] { energy_measurements_preserve_charged_state(scheme); });
            check::scenario(scheme + ": initial density uses the same charge shape as the step",
                            [&] { initial_density_uses_the_scheme_charge_shape(scheme); });
        }
        check::scenario("Corr velocity correction: electric work, invalid energy and unchanged failed state",
                        velocity_correction_adds_current_work_and_rejects_invalid_energy);
        check::scenario("ECSIM-Corr: energy and Gauss law across one/two periodic Z cells",
                        small_periodic_energy);
    });
}
