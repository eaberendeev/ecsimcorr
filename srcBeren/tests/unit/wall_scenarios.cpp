#include <filesystem>
#include <string>
#include <type_traits>

#include "simulation_ecsim_corr.h"
#include "test_support.h"

namespace {
constexpr int cells = 4;
constexpr double dt = 0.05;
constexpr double density = 1e-6;
constexpr double speed = 0.5;

nlohmann::json settings(const std::string& scheme, const nlohmann::json& boundaries) {
    return {{"Scheme", scheme},
            {"Dt", dt},
            {"NumCellsX", cells},
            {"NumCellsY", cells},
            {"NumCellsZ", cells},
            {"Dx", 1},
            {"Dy", 1},
            {"Dz", 1},
            {"StartFromTime", 0},
            {"StartTimeStep", 0},
            {"LastTimestep", 1},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collider", "None"},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions", boundaries}};
}

// One electron with a tiny macroparticle weight keeps the self-consistent
// fields negligible, so the wall test can compare against exact kinematics.
template <typename Scheme>
class WallScenario : public Scheme {
   public:
    using Scheme::Scheme;
    enum class WallCase { Absorbed, Reflected };

    WallCase wall_case = WallCase::Absorbed;
    Face wall = Face::XMAX;
    int samples = 0;

    void make_diagnostic(int step) override {
        Scheme::make_diagnostic(step);
        auto& electrons = *this->species.at("Electrons");
        const int count = electrons.get_total_num_of_particles();
        // Surviving macroparticle charge q*mpw*N. This is not the deposited
        // grid charge; the deposit is checked separately below and in
        // simulation_scenarios.
        const double particle_charge = electrons.charge * electrons.mpw() * count;

        if (step == 0) {
            check::require(count == 1, "wall scenario needs one initial electron");
            check::near(particle_charge, -density, 1e-18, "initial electron macroparticle charge");
            ++samples;
            return;
        }

        check::require(step == 1, "wall scenario runs exactly one step");
        if (wall_case == WallCase::Absorbed) {
            // The open wall removes the particle and accounts for its kinetic
            // energy, leaving no macroparticle charge behind.
            check::require(count == 0, "open wall must absorb the electron");
            const auto& lost = electrons.diag.boundary[wall];
            check::require(lost.lost_count == 1, "open wall must record one lost particle");
            check::require(lost.reflected_count == 0, "open wall must not reflect");
            check::near(lost.lost_energy, 0.5 * electrons.mpw() * electrons.mass() * speed * speed, 1e-15,
                        "open wall lost kinetic energy");
            check::near(particle_charge, 0, 0, "absorbed electron leaves no macroparticle charge");

            // With no particles left the deposit must be empty everywhere.
            // Sum only the interior nodes so ghost copies are not counted twice.
            Field3d density_grid(this->domain.size(), 1);
            this->collect_charge_density(density_grid);
            const int first = this->domain.grid.ghost_cells();
            double deposited_charge = 0;
            for (int i = first; i < first + cells; ++i)
                for (int j = first; j < first + cells; ++j)
                    for (int k = first; k < first + cells; ++k) deposited_charge += density_grid(i, j, k, 0);
            check::near(deposited_charge, 0, 1e-15, "absorbed electron must leave an empty density deposit");
        } else {
            // Reflection reverses the normal velocity and mirrors the endpoint,
            // conserving the particle and its macroparticle charge. Deposit
            // correctness at a reflecting wall is covered by simulation_scenarios.
            check::require(count == 1, "reflecting wall must keep the electron");
            const auto& reflected = electrons.diag.boundary[wall];
            check::require(reflected.reflected_count == 1, "reflecting wall must record one reflection");
            check::require(reflected.lost_count == 0, "reflecting wall must not lose the electron");
            bool found = false;
            for (int cell = 0; cell < electrons.size(); ++cell)
                for (const auto& particle : electrons.particlesData(cell)) {
                    found = true;
                    check::require(particle.coord.z() < cells, "reflected electron must stay inside the domain");
                    check::require(particle.velocity.z() < 0, "reflected electron must move back into the domain");
                    check::near(particle.velocity.z(), -speed, 1e-3, "reflected normal velocity");
                    check::near(particle.coord.z(), cells - 0.015, 1e-3, "reflected position after one step");
                }
            check::require(found, "reflected electron must remain on the grid");
            check::near(particle_charge, -density, 1e-18, "reflection conserves macroparticle charge");
        }
        ++samples;
    }
};

template <typename Scheme>
void exercise(const std::string& scheme, typename WallScenario<Scheme>::WallCase wall_case) {
    using Case = typename WallScenario<Scheme>::WallCase;
    const bool absorbed = wall_case == Case::Absorbed;
    // Open X walls absorb; the reflecting wall is a Z face because electron
    // reflection only supports ZMIN/ZMAX and the species named "Electrons".
    const char* absorbed_boundaries = R"([{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}},)"
                                      R"({"open":[{"face":"XMIN"},{"face":"XMAX"}]}])";
    const char* reflected_boundaries =
        R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"open":[{"face":"ZMIN"}]},)"
        R"({"electron_reflection":{"face":"ZMAX","radius":10,"energy_threshold":100}}])";
    const nlohmann::json boundaries = nlohmann::json::parse(absorbed ? absorbed_boundaries : reflected_boundaries);

    auto particles =
        nlohmann::json{{"particles", nlohmann::json::array({check::species("Electrons", -1, 1, density, 1)})}};
    char program[] = "wall_scenarios";
    char* argv[] = {program};
    WallScenario<Scheme> simulation(settings(scheme, boundaries), particles, 1, argv);
    simulation.wall_case = wall_case;
    // Electron reflection supports only Z faces; absorption uses the X wall.
    simulation.wall = absorbed ? Face::XMAX : Face::ZMAX;
    simulation.init();

    auto& electrons = *simulation.species.at("Electrons");
    const Vector3R position = absorbed ? Vector3R(cells - 0.01, 2, 2) : Vector3R(2, 2, cells - 0.01);
    const Vector3R velocity = absorbed ? Vector3R(speed, 0, 0) : Vector3R(0, 0, speed);
    electrons.add_particle(Particle(position, velocity));
    simulation.deposit_density(electrons);
    electrons.kineticEnergy = electrons.get_kinetic_energy();

    simulation.calculate();
    simulation.finalize();
    check::require(simulation.samples == 2, "wall scenario must be observed before and after the step");
}
}   // namespace

int main() {
    return check::run([] {
        const auto root = std::filesystem::current_path();
        using Case = WallScenario<SimulationEcsim>::WallCase;
        for (const std::string scheme : {"ecsim", "ecsim_corr"}) {
            for (Case wall_case : {Case::Absorbed, Case::Reflected}) {
                const bool absorbed = wall_case == Case::Absorbed;
                const std::string name = scheme + (absorbed ? "_open" : "_reflected");
                const auto directory = root / name;
                std::filesystem::create_directories(directory);
                std::filesystem::current_path(directory);
                if (scheme == "ecsim") {
                    check::scenario(name + ": full step at a wall",
                                    [&] { exercise<SimulationEcsim>(scheme, wall_case); });
                } else {
                    using CorrCase = WallScenario<SimulationEcsimCorr>::WallCase;
                    const CorrCase corr_case = absorbed ? CorrCase::Absorbed : CorrCase::Reflected;
                    check::scenario(name + ": full step at a wall",
                                    [&] { exercise<SimulationEcsimCorr>(scheme, corr_case); });
                }
                std::filesystem::current_path(root);
            }
        }
    });
}
