#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "simulation_ecsim.h"
#include "simulation_ecsim_corr.h"
#include "test_support.h"

namespace {
constexpr int cells = 4;
constexpr int split_step = 2;
constexpr int final_step = 4;
constexpr double dt = 0.2;
constexpr double tolerance = 1e-12;
using ParticleState = std::array<double, 6>;

struct EndpointState {
    Field3d electric;
    Field3d magnetic;
    std::map<std::string, std::vector<ParticleState>> particles;
};

nlohmann::json settings(const std::string& scheme, int first_step, int last_step, bool restart) {
    return {{"Scheme", scheme},
            {"Dt", dt},
            {"NumCellsX", cells},
            {"NumCellsY", cells},
            {"NumCellsZ", cells},
            {"Dx", 1},
            {"Dy", 1},
            {"Dz", 1},
            {"StartFromTime", restart ? 1 : 0},
            {"StartTimeStep", first_step},
            {"LastTimestep", last_step},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array({{{"type", "recovery"}, {"interval", 2}}})}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

template <typename Scheme>
class RestartPlasma : public Scheme {
   public:
    using Scheme::Scheme;
    std::map<int, EndpointState> endpoints;

    void make_diagnostic(int step) override {
        Scheme::make_diagnostic(step);
        EndpointState state{this->fieldEn, this->fieldBn, {}};
        for (const auto& [name, species] : this->species) {
            auto& particles = state.particles[name];
            for (int cell = 0; cell < species->size(); ++cell)
                for (const Particle& particle : (*species)(cell)) {
                    ParticleState values;
                    for (int axis = 0; axis < 3; ++axis) {
                        values[axis] = particle.coord[axis];
                        values[axis + 3] = particle.velocity[axis];
                        check::require(std::isfinite(values[axis]) && std::isfinite(values[axis + 3]),
                                       "restart plasma must have finite particle coordinates and velocities");
                    }
                    particles.push_back(values);
                }
            // Redistribution may change insertion order without changing the state.
            std::sort(particles.begin(), particles.end());
        }
        endpoints.emplace(step, std::move(state));
    }
};

template <typename Scheme>
void populate(RestartPlasma<Scheme>& simulation) {
    for (int i = 0; i < cells; ++i)
        for (int j = 0; j < cells; ++j)
            for (int k = 0; k < cells; ++k) {
                const Vector3R position(i + 0.99, j + 0.35, k + 0.65);
                const double phase = i + 2 * j + 3 * k + 0.5;
                const Vector3R velocity(0.04 * std::sin(phase), 0.04 * std::cos(phase), 0.03);
                for (auto& [name, species] : simulation.species)
                    for (int sign : {-1, 1})
                        species->add_particle(Particle(position, velocity * (sign / std::sqrt(species->mass()))));
            }
    for (auto& [name, species] : simulation.species)
        species->step_state.kinetic_energy_before_step = species->get_kinetic_energy();
}

void compare_endpoints(const EndpointState& actual, const EndpointState& expected, const std::string& label) {
    for (const auto pair :
         {std::pair{&actual.electric, &expected.electric}, std::pair{&actual.magnetic, &expected.magnetic}}) {
        check::equal(pair.first->capacity(), pair.second->capacity(), label + ": field size");
        for (int i = 0; i < pair.first->capacity(); ++i)
            check::near((*pair.first)[i], (*pair.second)[i], tolerance, label + ": endpoint field");
    }
    check::equal(actual.particles.size(), expected.particles.size(), label + ": species count");
    for (const auto& [name, particles] : expected.particles) {
        const auto& observed = actual.particles.at(name);
        check::equal(observed.size(), particles.size(), label + ": " + name + " count");
        for (size_t i = 0; i < particles.size(); ++i)
            for (int component = 0; component < 6; ++component)
                check::near(observed[i][component], particles[i][component], tolerance,
                            label + ": " + name + " coordinate/velocity");
    }
}

// The existing checkpoint stores endpoint E/B and particles. This experiment
// deliberately excludes processes whose RNG/accumulators are not saved:
// collisions, sources, emission and diagnostic averaging windows.
template <typename Scheme>
void continuous_matches_restarted(const std::string& scheme) {
    const auto root = std::filesystem::current_path();
    const auto directory = root / scheme;
    const auto continuous_directory = directory / "continuous";
    const auto split_directory = directory / "split";
    const auto resumed_directory = split_directory / "resumed";
    for (const auto& path : {continuous_directory, resumed_directory}) std::filesystem::create_directories(path);
    const nlohmann::json particles = {{"particles", nlohmann::json::array({check::species("Electrons", -1, 1, 1, 2),
                                                                           check::species("Ions", 1, 1836, 1, 2)})}};
    char program[] = "simulation_restart";
    char* argv[] = {program};

    std::filesystem::current_path(continuous_directory);
    RestartPlasma<Scheme> continuous(settings(scheme, 0, final_step, false), particles, 1, argv);
    continuous.init();
    populate(continuous);
    continuous.calculate();
    continuous.finalize();
    check::require(continuous.endpoints.at(final_step).electric.norm() > 1e-8,
                   "restart plasma must exercise nonzero self-consistent electric fields");

    std::filesystem::current_path(split_directory);
    RestartPlasma<Scheme> first(settings(scheme, 0, split_step, false), particles, 1, argv);
    first.init();
    populate(first);
    first.calculate();
    first.finalize();
    compare_endpoints(first.endpoints.at(split_step), continuous.endpoints.at(split_step),
                      scheme + ": independently initialized split prefix");

    // Production recovery reads ../Recovery. A new simulation reconstructs all
    // field operators and scratch storage before continuing the actual loop.
    std::filesystem::current_path(resumed_directory);
    RestartPlasma<Scheme> resumed(settings(scheme, split_step, final_step, true), particles, 1, argv);
    resumed.init();
    resumed.calculate();
    resumed.finalize();
    check::equal(resumed.endpoints.size(), size_t{1 + final_step - split_step},
                 "restart must observe initialization and every remaining step");
    compare_endpoints(resumed.endpoints.at(0), continuous.endpoints.at(split_step),
                      scheme + ": restored endpoint before continuation");
    for (int step = split_step + 1; step <= final_step; ++step)
        compare_endpoints(resumed.endpoints.at(step), continuous.endpoints.at(step),
                          scheme + ": resumed production step " + std::to_string(step));
    std::filesystem::current_path(root);
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("ECSIM: four continuous steps match two steps plus checkpoint/restart",
                        [] { continuous_matches_restarted<SimulationEcsim>("ecsim"); });
        check::scenario("ECSIM-Corr: four continuous steps match two steps plus checkpoint/restart",
                        [] { continuous_matches_restarted<SimulationEcsimCorr>("ecsim_corr"); });
    });
}
