#include <cmath>
#include <utility>
#include <vector>

#include "collisions/collision.h"
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

std::vector<Vector3R> velocities(const Species& species) {
    std::vector<Vector3R> result;
    for (const std::string name : {"Electrons", "Ions"}) {
        const auto& array = *species.at(name);
        for (int cell = 0; cell < array.size(); ++cell)
            for (const auto& particle : array(cell)) result.push_back(particle.velocity);
    }
    return result;
}

void same_velocities(const Species& actual, const std::vector<Vector3R>& expected, const std::string& message) {
    const auto observed = velocities(actual);
    check::equal(observed.size(), expected.size(), message + ": particle count");
    for (size_t i = 0; i < observed.size(); ++i)
        check::vector(observed[i], expected[i], 0, message + ": particle " + std::to_string(i));
}

bool velocities_differ(const Species& first, const Species& second) {
    const auto a = velocities(first), b = velocities(second);
    for (size_t i = 0; i < a.size(); ++i)
        if ((a[i] - b[i]).norm() > 1e-14)
            return true;
    return false;
}

struct Moments {
    Vector3R momentum{0, 0, 0};
    double energy = 0;
};

Moments moments(const Species& species) {
    Moments result;
    for (const auto& entry : species) {
        const auto& array = *entry.second;
        for (int cell = 0; cell < array.size(); ++cell)
            for (const auto& particle : array(cell)) {
                const double mass = array.mpw() * array.mass();
                result.momentum += mass * particle.velocity;
                result.energy += 0.5 * mass * particle.velocity.dot(particle.velocity);
            }
    }
    return result;
}

void conserved(const Species& species, const Moments& initial, const std::string& message) {
    const auto final = moments(species);
    check::near(final.energy, initial.energy, 1e-14, message + ": kinetic energy");
    check::vector(final.momentum, initial.momentum, 1e-13, message + ": momentum");
}

nlohmann::json coulomb_channel(const std::string& first, const std::string& second) {
    return {{"type", "coulomb"}, {"species1", first}, {"species2", second}};
}

void disabled_and_rejected_configuration(const Domain& domain) {
    auto species = mixture(domain);
    const auto initial = velocities(species);
    auto disabled = collisions_from_config({{"Collider", "None"}}, density_normalization);
    check::require(!disabled.has_collisions() && !disabled.enabled(),
                   "missing modern channels and None disable collisions");
    check::equal(disabled.legacy_info(), std::string("None"), "disabled legacy description");
    disabled.apply(species, domain, dt);
    same_velocities(species, initial, "disabled manager preserves velocities");

    // Missing/wrong legacy settings retain their explicit configuration error when no channel was accepted.
    for (const auto& config : {nlohmann::json::object(), nlohmann::json{{"Collider", 3}}}) {
        bool rejected = false;
        try {
            auto manager = collisions_from_config(config, density_normalization);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check::require(rejected, "missing or non-string Collider must be rejected without modern channels");
    }
    for (const auto& entry :
         {nlohmann::json{{"type", 3}}, nlohmann::json{{"type", "coulomb"}, {"species1", 3}, {"species2", "Ions"}}}) {
        bool rejected = false;
        try {
            auto manager = collisions_from_config(
                {{"Collider", "None"}, {"Collisions", nlohmann::json::array({entry})}}, density_normalization);
        } catch (const nlohmann::json::exception&) {
            rejected = true;
        }
        check::require(rejected, "incorrect modern field types retain JSON type errors");
    }
}

void unaccepted_channels_fall_back(const Domain& domain) {
    const std::vector<nlohmann::json> unaccepted{
        nlohmann::json::array(), 7,
        nlohmann::json::array({nlohmann::json::object(),
                               {{"type", "unknown"}},
                               {{"type", "coulomb"}, {"species1", "Electrons"}},
                               {{"type", "neutral_ionization"}}})};
    for (const auto& channels : unaccepted) {
        auto manager =
            collisions_from_config({{"Collider", "Coulomb"}, {"Collisions", channels}}, density_normalization);
        check::require(!manager.has_collisions() && manager.enabled(), "unaccepted channels select legacy fallback");
        check::equal(manager.legacy_info(), std::string("Coulomb"), "legacy fallback name");
        auto actual = mixture(domain), expected = mixture(domain);
        BinaryCollider reference(density_normalization);
        manager.apply(actual, domain, dt);
        reference.collide_ion_electron_binary(expected, dt);
        same_velocities(actual, velocities(expected), "unaccepted modern channels execute legacy Coulomb");
    }
    auto unknown = collisions_from_config({{"Collider", "UnknownLegacy"}}, density_normalization);
    auto species = mixture(domain);
    const auto initial = velocities(species);
    unknown.apply(species, domain, dt);
    same_velocities(species, initial, "unknown legacy name retains no-op behavior");
    check::equal(unknown.legacy_info(), std::string("UnknownLegacy"),
                 "unknown legacy name remains available for logging");
}

void modern_order_defaults_and_move(const Domain& domain) {
    const auto same = coulomb_channel("Electrons", "Electrons");
    auto different = coulomb_channel("Electrons", "Ions");
    different["coulomb_log"] = 7.0;
    const nlohmann::json channels = nlohmann::json::array({same, different});
    // A valid modern list works without Collider, and also takes precedence over legacy Coulomb.
    auto manager = collisions_from_config({{"Collisions", channels}}, density_normalization);
    auto independent =
        collisions_from_config({{"Collisions", channels}, {"Collider", "Coulomb"}}, density_normalization);
    check::require(manager.has_collisions() && manager.enabled(), "accepted modern channels enable collisions");
    check::equal(independent.legacy_info(), std::string("None"), "modern channels suppress legacy fallback");
    auto actual = mixture(domain), expected = mixture(domain), delayed = mixture(domain), reversed = mixture(domain);
    const auto initial_velocities = velocities(actual);
    const auto initial = moments(actual);
    same_velocities(actual, initial_velocities, "factory construction preserves initial velocities");
    CoulombCollisionOperator same_reference("Electrons", "Electrons", density_normalization, 15.0);
    CoulombCollisionOperator different_reference("Electrons", "Ions", density_normalization, 7.0);
    CoulombCollisionOperator reverse_same("Electrons", "Electrons", density_normalization, 15.0);
    CoulombCollisionOperator reverse_different("Electrons", "Ions", density_normalization, 7.0);
    for (int step = 0; step < 2; ++step) {
        manager.apply(actual, domain, dt);
        same_reference.apply(expected, domain, dt);
        different_reference.apply(expected, domain, dt);
        same_velocities(actual, velocities(expected), "configured order, default lnL=15 and per-channel streams");
        reverse_different.apply(reversed, domain, dt);
        reverse_same.apply(reversed, domain, dt);
    }
    check::require(velocities_differ(actual, reversed), "fixture must distinguish the two physical channel orders");
    auto moved = std::move(manager);
    for (int step = 2; step < 5; ++step) {
        moved.apply(actual, domain, dt);
        same_reference.apply(expected, domain, dt);
        different_reference.apply(expected, domain, dt);
        same_velocities(actual, velocities(expected), "move preserves each channel's apply index between steps");
        conserved(actual, initial, "ordered modern channels conserve");
    }
    // Applying another manager after the first has advanced must still start at that manager's own step zero.
    for (int step = 0; step < 5; ++step) independent.apply(delayed, domain, dt);
    same_velocities(delayed, velocities(actual), "modern managers own independent random streams");
    auto untouched = mixture(domain);
    check::require(velocities_differ(actual, untouched), "positive modern collision rate must scatter particles");
}

void disabled_neutral_reactions_preserve_particles_and_drift(const Domain& domain) {
    const nlohmann::json channel = {{"type", "neutral_ionization"},     {"charged_species", "Projectile"},
                                    {"neutral_species", "HydrogenGas"}, {"electron_product", "SecondaryElectrons"},
                                    {"ion_product", "ProductIons"},     {"scheme", "physical_only"},
                                    {"electron_ionization", false},     {"proton_ionization", false},
                                    {"proton_charge_exchange", false}};
    auto manager = collisions_from_config({{"Collisions", nlohmann::json::array({channel})}}, density_normalization);
    check::require(manager.has_collisions(), "custom neutral channel is accepted without legacy Collider");
    Species species;
    species.emplace("Projectile", std::make_unique<ParticlesArray>(check::species("Projectile", 1, 1836), domain));
    species.emplace("HydrogenGas", std::make_unique<ParticlesArray>(check::species("HydrogenGas", 0, 1836), domain));
    species.emplace("SecondaryElectrons",
                    std::make_unique<ParticlesArray>(check::species("SecondaryElectrons", -1, 1), domain));
    species.emplace("ProductIons", std::make_unique<ParticlesArray>(check::species("ProductIons", 1, 1836), domain));
    const Vector3R position(0.5, 0.5, 0.5), projectile_velocity(0.02, -0.01, 0.03);
    const Vector3R neutral_velocity(0.125, -0.0625, 0.03125);
    species.at("Projectile")->add_particle(Particle(position, projectile_velocity));
    species.at("HydrogenGas")->add_particle(Particle(position, neutral_velocity));
    // Dyadic dt and velocities make the two free drifts exactly representable.
    constexpr double neutral_dt = 0.125;
    for (int step = 0; step < 2; ++step) manager.apply(species, domain, neutral_dt);
    for (const auto& name : {"Projectile", "HydrogenGas"})
        check::equal(species.at(name)->get_total_num_of_particles(), 1,
                     std::string(name) + " survives disabled reactions");
    for (const auto& name : {"SecondaryElectrons", "ProductIons"})
        check::equal(species.at(name)->get_total_num_of_particles(), 0, std::string(name) + " receives no products");
    for (const auto& entry : species) {
        const auto& array = *entry.second;
        for (int cell = 0; cell < array.size(); ++cell)
            for (const auto& particle : array(cell)) {
                const bool neutral = entry.first == "HydrogenGas";
                check::vector(particle.velocity, neutral ? neutral_velocity : projectile_velocity, 0,
                              entry.first + " velocity is preserved with zero reaction rate");
                check::vector(particle.coord, neutral ? position + 2 * neutral_dt * neutral_velocity : position, 0,
                              entry.first + " has the expected two-step drift");
            }
    }
}

void legacy_steps_move_and_independence(const Domain& domain) {
    auto manager = collisions_from_config({{"Collider", "Coulomb"}}, density_normalization);
    auto independent = collisions_from_config({{"Collider", "Coulomb"}}, density_normalization);
    auto actual = mixture(domain), expected = mixture(domain), delayed = mixture(domain);
    const auto initial = moments(actual);
    const auto initial_velocities = velocities(expected);
    same_velocities(actual, initial_velocities, "legacy manager construction preserves initial velocities");
    BinaryCollider reference(density_normalization);
    manager.apply(actual, domain, dt);
    reference.collide_ion_electron_binary(expected, dt);
    same_velocities(actual, velocities(expected), "legacy first step retains BinaryCollider stream");
    auto moved = std::move(manager);
    for (int step = 1; step < 5; ++step) {
        moved.apply(actual, domain, dt);
        reference.collide_ion_electron_binary(expected, dt);
        same_velocities(actual, velocities(expected), "legacy repeated steps and move retain BinaryCollider stream");
        conserved(actual, initial, "legacy electron-ion collisions conserve");
    }
    for (int step = 0; step < 5; ++step) independent.apply(delayed, domain, dt);
    same_velocities(delayed, velocities(actual), "legacy managers own independent engines");
    auto untouched = mixture(domain);
    check::require(velocities_differ(actual, untouched), "positive legacy collision rate must scatter particles");
}
}   // namespace

int main() {
    return check::run([] {
        Domain domain;
        domain.init(Vector3I(4, 2, 2), Vector3R(1, 1, 1));
        check::scenario("collision configuration: disabled and rejected settings",
                        [&] { disabled_and_rejected_configuration(domain); });
        check::scenario("collision configuration: unaccepted modern channels use legacy fallback",
                        [&] { unaccepted_channels_fall_back(domain); });
        check::scenario("collision configuration: modern order, defaults and persistent channel streams",
                        [&] { modern_order_defaults_and_move(domain); });
        check::scenario("collision configuration: disabled custom neutral reactions retain free drift",
                        [&] { disabled_neutral_reactions_preserve_particles_and_drift(domain); });
        check::scenario("collision configuration: legacy steps, move and independent engines",
                        [&] { legacy_steps_move_and_independence(domain); });
    });
}
