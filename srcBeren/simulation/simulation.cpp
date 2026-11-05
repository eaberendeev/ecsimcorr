// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#include "Coil.h"
#include "Damping.h"
#include "Diagnostic.h"
#include "config/boundary_config.h"
#include "config/collision_config.h"
#include "config/domain_config.h"
#include "config/json_settings.h"
#include "config/particle_config.h"
#include "containers.h"
#include "grid/domain.h"
#include "io/field_io.h"
#include "io/recovery.h"
#include "io/text_input.h"
#include "particle_array.h"
#include "runtime/env_options.h"
#include "runtime/log_macros.h"
#include "runtime/profiler.h"
#include "simulation_ecsim.h"
#include "simulation_ecsim_corr.h"
#include "util.h"

Simulation::Simulation(const nlohmann::json& s_config, const nlohmann::json& p_config, int argc, char** argv)
    : system_config(s_config), particles_config(p_config) {
    (void) argv[argc - 1];
}

Simulation::Simulation() = default;

Simulation::~Simulation() = default;

void Simulation::init() {
    timer::profiling_enabled = true;
    RECORD_TIMER;
    dt_ = get_checked<double>(system_config, "Dt");
    domain = domain_from_config(system_config);
    bc_handler = boundary_conditions_from_config(system_config, domain);

    init_fields();
    init_particles(particles_config);
    // Проверка сортов-продуктов secondary_emission после того, как созданы все
    // сорта частиц (опечатка в "product" должна упасть с понятным сообщением
    // до основного цикла; здесь последовательный код — бросать можно).
    bc_handler.validate_emissions(species);

    const double n0 = get_checked<double>(system_config, "n0");
    collision_manager = collisions_from_config(system_config, n0);

    logger::init();
    g_verbose_step = system_config.value("VerboseStep", false);
    if (g_verbose_step)
        logger::info("VerboseStep enabled");

    std::cout << "Simulation initialized\n";
}

void Simulation::calculate() {
    const int startTimeStep = get_checked<int>(system_config, "StartTimeStep");
    const int lastTimestep = get_checked<int>(system_config, "LastTimestep");

    std::cout << "\n=== beren3d ===\n";
    std::cout << "Scheme: " << get_checked<std::string>(system_config, "Scheme") << "\n";
    std::cout << "Grid:  " << get_checked<int>(system_config, "NumCellsX") << "x"
              << get_checked<int>(system_config, "NumCellsY") << "x" << get_checked<int>(system_config, "NumCellsZ")
              << "  dx=" << get_checked<double>(system_config, "Dx") << "  dt=" << dt_ << "\n";
    std::cout << "Steps: " << startTimeStep + 1 << ".." << lastTimestep << "\n";

    if (collision_manager.enabled())
        std::cout << "Collisions: enabled\n";
    std::cout << "\n";

    const int timeredSimSteps = envOptions::timeredSimSteps();
    const int timerStartStep = envOptions::timerStartStep();
    const int timerEndStep = envOptions::timerEndStep();
    timer::profiling_enabled = timerStartStep <= 0 && timerEndStep >= 0;
    if (!timer::profiling_enabled)
        timer::clearFullProfile();
    make_diagnostic(0);
    for (auto timestep = startTimeStep + 1; timestep <= lastTimestep; ++timestep) {
        const bool inTimerWindow = timestep >= timerStartStep && timestep <= timerEndStep;
        timer::profiling_enabled = inTimerWindow;
        timer::clearTree();
        // Production wall time includes preparation, physics and diagnostic I/O.
        // Detailed profile serialization and the timing log message are excluded.
        const auto step_start = std::chrono::steady_clock::now();
        logger::set_timestep(timestep);
        timer::timer timerSimLoop("sim loop");

        timer::timer timerPrepare("Prepare step");
        prepare_step(timestep);
        make_step(timestep);
        timerPrepare.finish();

        timer::timer timerCollision("collision step");
        if (collision_manager.enabled()) {
            collision_step(timestep);
        }
        timerCollision.finish();

        timer::timer timerDiagnostic("diagnostic");
        make_diagnostic(timestep);
        timerDiagnostic.finish();
        const double wall_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - step_start).count();
        LOG_STEP("Step wall time: " << std::setprecision(17) << wall_seconds << " s");

        timerSimLoop.finish();
        if (!inTimerWindow) {
            // Disabled timer constructors skip clocks and event collection outside the window.
            timer::clearFullProfile();
        } else if (timeredSimSteps == 0) {
            timer::commonTimer timerWriteProfile("writing profile");
            timer::writeFullProfile("profile.json");
            timerWriteProfile.finish();
        } else if (timestep == timerEndStep || timestep % timeredSimSteps == 0) {
            timer::commonTimer timerWriteProfile("writing profile");
            timer::writeFullProfile("profile_it_" + std::to_string(timestep) + ".json");
            timerWriteProfile.finish();
            timer::clearFullProfile();
        }
    }
}

void Simulation::finalize() {
    logger::close();
}

void Simulation::init_particles(const nlohmann::json& j) {
    RECORD_TIMER;
    if (!j.contains("particles") || !j["particles"].is_array()) {
        throw std::runtime_error("mixture requires particles[]");
    }

    for (const nlohmann::json& config : j["particles"]) {
        std::unique_ptr<ParticlesArray> sp = make_particles_array(config);
        ParticleSources sources = particle_sources_from_config(config, sp->parameters(), domain);
        const std::string name = sp->name();
        if (name.empty()) {
            throw std::runtime_error("Particle species must have non-empty Name");
        }
        if (species.count(name)) {
            throw std::runtime_error("Duplicate particle species name: " + name);
        }
        particle_sources.emplace(name, std::move(sources));
        species.emplace(name, std::move(sp));
    }

    charged_species.clear();
    for (auto& kv : species) {
        ParticlesArray& sp = *kv.second;
        if (!sp.is_neutral()) {
            charged_species.emplace(kv.first, std::ref(sp));
        }
    }

    for (auto& kv : species) {
        ParticlesArray& sp = *kv.second;
        if (get_checked<double>(system_config, "k_particles_reservation") > 0.) {
            for (auto k = 0; k < sp.size(); ++k) {
                sp(k).reserve(get_checked<double>(system_config, "k_particles_reservation") * sp.particles_per_cell());
            }
        }
        if (get_checked<int>(system_config, "StartFromTime") > 0) {
            read_particles_from_recovery(&sp);
            std::cout << "  Loaded " << sp.name() << " from recovery\n";
        } else {
            double init_energy = distribute_initial_particles(sp, particle_sources.at(sp.name()).initial, domain,
                                                              system_config.value("co_locate_species", true));
            std::cout << "  " << sp.name() << " distributed, E_kin=" << init_energy << "\n";
        }
        sp.step_state.kinetic_energy_before_step = sp.get_kinetic_energy();
    }
}

void Simulation::advance_neutral_particles(double dt) {
    for (auto& entry : species) {
        auto& sp = *entry.second;
        if (!sp.is_neutral())
            continue;
        sp.move(dt);
        bc_handler.apply_to_particles(sp, species, domain);
        sp.update_cells(domain);
    }
}

void Simulation::collect_current(Field3d& J) {
    RECORD_TIMER;
    J.setZero();
    for (const auto& kv : species) {
        const auto& sp = *kv.second;
        if (!sp.is_neutral())
            J.data() += sp.currentOnGrid.data();
    }
}

std::unique_ptr<Simulation> build_simulation(const nlohmann::json& system_config,
                                             const nlohmann::json& particles_config, int argc, char** argv) {
    RECORD_TIMER;
    auto scheme_name = get_checked<std::string>(system_config, "Scheme");

    std::unique_ptr<Simulation> simulation = nullptr;

    if (scheme_name == "ecsim") {
        simulation = std::make_unique<SimulationEcsim>(system_config, particles_config, argc, argv);
    } else if (scheme_name == "ecsim_corr") {
        simulation = std::make_unique<SimulationEcsimCorr>(system_config, particles_config, argc, argv);
    } else {
        std::cout << "Scheme " << scheme_name << " is not supported\n";
        exit(-1);
    }
    return simulation;
}

void Simulation::collision_step([[maybe_unused]] const int timestep) {
    RECORD_TIMER;

    collision_manager.apply(species, dt_);
}
