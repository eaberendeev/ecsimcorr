#include <filesystem>
#include <fstream>
#include <iostream>

#include "simulation/eccapfim/algorithms_eccapfim.h"
#include "io/field_io.h"
#include "simulation/simulation.h"

// Replay a single trial orbit with the exact total Lorentz fields from a
// failure dump; no injection, collisions or nonlinear field solve is needed.
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: eccapfim_replay <EccapfimFailure/step_N>\n";
        return 2;
    }
    try {
        const auto directory = std::filesystem::absolute(argv[1]);
        std::ifstream file(directory / "failure.json");
        nlohmann::json report;
        file >> report;
        auto config = report.at("system_config");
        auto particles = report.at("particles_config");
        config["StartFromTime"] = 0;
        config["StartTimeStep"] = 0;
        config["Collisions"] = nlohmann::json::array();
        config["diagnostics"]["outputs"] = nlohmann::json::array();
        config["Eccapfim"]["failure_dump"] = false;
        config["Eccapfim"]["iteration_log"] = false;
        for (auto& species : particles.at("particles")) species["distribution"] = nlohmann::json::array();
        char program[] = "eccapfim_replay";
        char* arguments[] = {program};
        auto simulation = build_simulation(config, particles, 1, arguments);
        simulation->init();
        const auto& domain = simulation->domain;
        const auto& bc = simulation->bc_handler;
        const auto& species = *simulation->species.at(report.at("species").get<std::string>());
        Field3d electric(domain.size(), 3), magnetic(domain.size(), 3), current(domain.size(), 3);
        read_field_from_file((directory / "electric.backup").string(), electric);
        read_field_from_file((directory / "magnetic.backup").string(), magnetic);
        current.setZero();
        const auto vector = [](const nlohmann::json& v) { return Vector3R(v.at(0), v.at(1), v.at(2)); };
        const auto json_vector = [](const Vector3R& v) { return nlohmann::json::array({v.x(), v.y(), v.z()}); };
        const auto& solver = report.at("particle_solver");
        const Particle initial(vector(solver.at("initial").at("position")),
                               vector(solver.at("initial").at("velocity")));
        const double dt = solver.at("dt");
        algorithmsECCAPFIM::OrbitBoundaryEvents events;
        const auto end = algorithmsECCAPFIM::advance_particle(
            initial, electric, magnetic, current, solver.at("charge_mass"), solver.at("charge_weight"), dt, domain,
            solver.at("tolerance"), solver.at("max_iterations"), &bc, bc.path_policy(species, domain), &events);
        // Count periodic owner nodes once; open axes retain all storage nodes.
        IndexRange owners;
        for (int d = 0; d < 3; ++d) {
            owners.start[d] = bc.is_periodic(d) ? GHOST_CELLS : 0;
            owners.end[d] = bc.is_periodic(d) ? owners.start[d] + domain.num_cells(d) : domain.size()[d];
        }
        const double work = dt * dot_product_sum(electric, current, owners);
        const double kinetic_change =
            0.5 * species.mass() * species.mpw() * (end.velocity.squaredNorm() - initial.velocity.squaredNorm());
        const nlohmann::json result{{"source_step", report.at("step")},
                                    {"position", json_vector(end.coord)},
                                    {"velocity", json_vector(end.velocity)},
                                    {"removed_face", static_cast<int>(events.removed)},
                                    {"reflections", events.reflected},
                                    {"kinetic_change", kinetic_change},
                                    {"electric_work", work},
                                    {"energy_residual", kinetic_change - work}};
        std::cout << result.dump(2) << std::endl;
        simulation->finalize();
        return 0;
    } catch (const algorithmsECCAPFIM::ParticleSolveFailure& error) {
        std::cerr << error.details.dump(2) << std::endl;
        return 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 2;
    }
}
