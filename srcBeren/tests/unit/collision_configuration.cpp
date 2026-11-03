#include <cmath>
#include <utility>
#include <vector>

#include "collisions/coulomb.h"
#include "config/collision_config.h"
#include "test_support.h"

namespace {
constexpr double density_normalization = 1e10;
constexpr double dt = 0.03;

Species mixture(const Domain& domain) {
    Species species;
    for (const std::string name : {"Electrons", "Ions"}) {
        const bool electron = name == "Electrons";
        auto array =
            std::make_unique<ParticlesArray>(check::species(name, electron ? -1 : 1, electron ? 1 : 7), domain);
        for (int cell = 0; cell < 3; ++cell)
            for (int i = 0; i < (electron ? 5 : 3); ++i) {
                const double phase = i + cell * 0.7 + (electron ? 0.1 : 0.6);
                array->add_particle(Particle(Vector3R(0.5 + cell, 0.5, 0.5),
                                             Vector3R(0.01 * std::sin(phase), 0.02 * std::cos(phase), 0.01 * (i - 2))));
            }
        species.emplace(name, std::move(array));
    }
    return species;
}

std::vector<Particle> particles(const Species& species) {
    std::vector<Particle> result;
    for (const std::string name : {"Electrons", "Ions"}) {
        const auto& array = *species.at(name);
        for (int cell = 0; cell < array.size(); ++cell)
            for (const auto& particle : array(cell)) result.push_back(particle);
    }
    return result;
}

void same_particles(const Species& actual, const std::vector<Particle>& expected, const std::string& message) {
    const auto observed = particles(actual);
    check::equal(observed.size(), expected.size(), message + ": particle count");
    for (size_t i = 0; i < observed.size(); ++i) {
        const auto label = message + ": particle " + std::to_string(i);
        check::vector(observed[i].velocity, expected[i].velocity, 0, label + " velocity");
        check::vector(observed[i].coord, expected[i].coord, 0, label + " coordinate");
    }
}

bool velocities_differ(const Species& first, const Species& second) {
    const auto a = particles(first), b = particles(second);
    for (size_t i = 0; i < a.size(); ++i)
        if ((a[i].velocity - b[i].velocity).norm() > 1e-14)
            return true;
    return false;
}

nlohmann::json coulomb_channel(const std::string& first, const std::string& second) {
    return {{"type", "coulomb"}, {"species1", first}, {"species2", second}};
}

void disabled_configuration(const Domain& domain) {
    const std::vector<nlohmann::json> configurations{
        nlohmann::json::object(),
        {{"Collisions", nlohmann::json::array()}},
        {{"Collisions", 7}},
        {{"Collisions", nlohmann::json::array({nlohmann::json::object(),
                                               {{"type", "unknown"}},
                                               {{"type", "coulomb"}, {"species1", "Electrons"}},
                                               {{"type", "neutral_ionization"}}})}}};
    for (const auto& config : configurations) {
        auto manager = collisions_from_config(config, density_normalization);
        check::require(!manager.enabled(), "missing, empty or skipped channels disable collisions");
        auto species = mixture(domain);
        const auto initial = particles(species);
        manager.apply(species, dt);
        same_particles(species, initial, "disabled manager preserves all particles");
    }
}

void rejected_configuration() {
    const auto valid_channels = nlohmann::json::array({coulomb_channel("Electrons", "Ions")});
    for (const auto& selector :
         {nlohmann::json("None"), nlohmann::json("Coulomb"), nlohmann::json(3), nlohmann::json(nullptr)}) {
        for (auto config : {nlohmann::json::object(), nlohmann::json{{"Collisions", nlohmann::json::array()}},
                            nlohmann::json{{"Collisions", valid_channels}}}) {
            config["Collider"] = selector;
            bool rejected = false;
            try {
                auto manager = collisions_from_config(config, density_normalization);
            } catch (const std::invalid_argument& error) {
                const std::string message = error.what();
                rejected =
                    message.find("Collider") != std::string::npos && message.find("Collisions") != std::string::npos;
            }
            check::require(rejected, "deprecated Collider requires a migration error even with modern channels");
        }
    }
    for (const auto& entry :
         {nlohmann::json{{"type", 3}}, nlohmann::json{{"type", "coulomb"}, {"species1", 3}, {"species2", "Ions"}},
          nlohmann::json{
              {"type", "coulomb"}, {"species1", "Electrons"}, {"species2", "Ions"}, {"coulomb_log", "invalid"}},
          nlohmann::json{
              {"type", "neutral_ionization"}, {"charged_species", "Projectile"}, {"electron_ionization", "invalid"}}}) {
        bool rejected = false;
        try {
            auto manager =
                collisions_from_config({{"Collisions", nlohmann::json::array({entry})}}, density_normalization);
        } catch (const nlohmann::json::exception&) {
            rejected = true;
        }
        check::require(rejected, "incorrect modern field types retain JSON type errors");
    }
}

void modern_order_defaults_and_streams(const Domain& domain) {
    const auto same = coulomb_channel("Electrons", "Electrons");
    auto different = coulomb_channel("Electrons", "Ions");
    different["coulomb_log"] = 7.0;
    // Skipped entries do not reorder accepted channels.
    const nlohmann::json channels = nlohmann::json::array({nlohmann::json::object(),
                                                           same,
                                                           {{"type", "unknown"}},
                                                           different,
                                                           {{"type", "coulomb"}, {"species1", "Electrons"}},
                                                           {{"type", "neutral_ionization"}}});
    auto manager = collisions_from_config({{"Collisions", channels}}, density_normalization);
    auto independent = collisions_from_config({{"Collisions", channels}}, density_normalization);
    check::require(manager.enabled(), "accepted modern channels enable collisions");
    auto actual = mixture(domain), expected = mixture(domain), delayed = mixture(domain), reversed = mixture(domain);
    CoulombCollisionOperator same_reference("Electrons", "Electrons", density_normalization, 15.0);
    CoulombCollisionOperator different_reference("Electrons", "Ions", density_normalization, 7.0);
    CoulombCollisionOperator reverse_same("Electrons", "Electrons", density_normalization, 15.0);
    CoulombCollisionOperator reverse_different("Electrons", "Ions", density_normalization, 7.0);
    for (int step = 0; step < 5; ++step) {
        manager.apply(actual, dt);
        same_reference.apply(expected, dt);
        different_reference.apply(expected, dt);
        same_particles(actual, particles(expected), "configured order, default lnL=15 and per-channel streams");
        reverse_different.apply(reversed, dt);
        reverse_same.apply(reversed, dt);
    }
    check::require(velocities_differ(actual, reversed), "fixture must distinguish the two physical channel orders");
    // Applying another manager after the first has advanced must still start at that manager's own step zero.
    for (int step = 0; step < 5; ++step) independent.apply(delayed, dt);
    same_particles(delayed, particles(actual), "modern managers own independent random streams");
    auto untouched = mixture(domain);
    check::require(velocities_differ(actual, untouched), "positive modern collision rate must scatter particles");
}

void neutral_channels_do_not_transport_particles(const Domain& domain) {
    const std::string projectile = "Projectile", gas = "HydrogenGas";
    const std::string electrons = "SecondaryElectrons", ions = "ProductIons";
    const nlohmann::json channel = {
        {"type", "neutral_ionization"},  {"charged_species", projectile}, {"neutral_species", gas},
        {"electron_product", electrons}, {"ion_product", ions},           {"scheme", "physical_only"},
        {"electron_ionization", false},  {"proton_ionization", false},    {"proton_charge_exchange", false}};
    auto manager =
        collisions_from_config({{"Collisions", nlohmann::json::array({channel, channel})}}, density_normalization);
    check::require(manager.enabled(), "two reaction-disabled custom neutral channels are enabled");
    Species species;
    for (const auto& name : {projectile, electrons, ions, gas}) {
        const bool neutral = name == gas;
        const bool electron = name == electrons;
        species.emplace(name,
                        std::make_unique<ParticlesArray>(
                            check::species(name, neutral ? 0 : (electron ? -1 : 1), electron ? 1 : 1836), domain));
    }
    const Vector3R position(0.5, 0.5, 0.5), gas_velocity(0.125, -0.0625, 0.03125);
    const Vector3R projectile_velocity(0.02, -0.01, 0.03);
    species.at(projectile)->add_particle(Particle(position, projectile_velocity));
    species.at(gas)->add_particle(Particle(position, gas_velocity));
    // With every reaction disabled, collision channels must not advance either species in space.
    for (int step = 0; step < 2; ++step) {
        manager.apply(species, 0.125);
        for (const auto& entry : species) {
            const bool populated = entry.first == projectile || entry.first == gas;
            const auto& array = *entry.second;
            const auto label = "reaction-disabled channels: " + entry.first + " step " + std::to_string(step);
            check::equal(array.get_total_num_of_particles(), populated ? 1 : 0, label + " particle count");
            for (int cell = 0; cell < array.size(); ++cell)
                for (const auto& particle : array(cell)) {
                    check::vector(particle.velocity, entry.first == gas ? gas_velocity : projectile_velocity, 0,
                                  label + " velocity");
                    check::vector(particle.coord, position, 0, label + " coordinate");
                }
        }
    }
}
}   // namespace

int main() {
    return check::run([] {
        Domain domain;
        domain.init(Vector3I(4, 2, 2), Vector3R(1, 1, 1));
        check::scenario("collision configuration: missing, empty and skipped channels disable collisions",
                        [&] { disabled_configuration(domain); });
        check::scenario("collision configuration: deprecated selector and wrong modern field types",
                        [] { rejected_configuration(); });
        check::scenario("collision configuration: modern order, defaults and persistent channel streams",
                        [&] { modern_order_defaults_and_streams(domain); });
        check::scenario("collision configuration: neutral channels do not transport particles",
                        [&] { neutral_channels_do_not_transport_particles(domain); });
    });
}
