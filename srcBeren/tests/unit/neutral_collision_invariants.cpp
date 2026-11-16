#include <cmath>
#include <iomanip>
#include <memory>
#include <vector>

#include "collisions/neutral/collisions_with_neutrals.h"
#include "collisions/neutral_collision.h"
#include "test_support.h"

void zero_limits_and_charge_exchange() {
    CollisionProcessOptions options;
    options.electron_ionization = false;
    options.proton_ionization = false;
    options.proton_charge_exchange = false;
    ColliderWithNeutrals disabled(1e10, CollisionScheme::PHYSICAL_ONLY, options);
    Vector3R cp(0.01, -0.002, 0.003), neutral(-0.001, 0.002, 0.001);
    check::near(disabled.total_collision_frequency(cp, neutral, 1836, 1), 0, 0, "disabled channels still scatter");
    for (CollisionScheme scheme : {CollisionScheme::PHYSICAL_ONLY, CollisionScheme::NULL_COLLISION}) {
        options.proton_charge_exchange = true;
        ColliderWithNeutrals collider(1e10, scheme, options);
        collider.rng.seed(17);
        check::near(collider.total_collision_frequency(cp, neutral, 1836, 0), 0, 0, "vacuum collision rate");
        auto a = cp, b = neutral;
        auto [ionized, va, vb] = collider.collision_with_neutral(a, b, 1836, 1836, 1, 0, 1);
        check::require(!ionized, "zero timestep created particles");
        check::vector(va, cp, 0, "zero timestep changed charged velocity");
        check::vector(vb, neutral, 0, "zero timestep changed neutral velocity");
        const double frequency = collider.total_collision_frequency(cp, neutral, 1836, 1e20);
        check::require(std::isfinite(frequency) && frequency > 0, "enabled CX channel must have a rate");
        a = cp;
        b = neutral;
        std::tie(ionized, va, vb) = collider.collision_with_neutral(a, b, 1836, 1836, 1e20, 1, frequency);
        check::require(!ionized, "charge exchange must not create an ionization pair");
        check::vector(va, neutral, 0, "charge exchange must swap charged/neutral velocities");
        check::vector(vb, cp, 0, "charge exchange must conserve pair momentum and energy");
    }
}

// A saturated proposal is accepted with probability one. With nu_max = 2*nu_CX,
// channel selection alone must produce half CX and half null events.
void seeded_null_channel_selection() {
    CollisionProcessOptions options;
    options.electron_ionization = false;
    options.proton_ionization = false;
    ColliderWithNeutrals collider(1e10, CollisionScheme::NULL_COLLISION, options);
    const Vector3R incident(0.01, 0, 0), stationary(0.0, 0.0, 0.0);
    const double frequency = collider.total_collision_frequency(incident, stationary, 1836, 1);
    check::require(std::isfinite(frequency) && frequency > 0, "CX rate for null-channel sampling");
    constexpr int trials = 4096;
    std::vector<bool> first_run;
    int reactions = 0;
    for (int replay = 0; replay < 2; ++replay) {
        collider.rng.seed(1729);
        for (int trial = 0; trial < trials; ++trial) {
            Vector3R projectile = incident, target = stationary;
            const auto [ionized, va, vb] =
                collider.collision_with_neutral(projectile, target, 1836, 1836, 1, 100 / frequency, 2 * frequency);
            check::require(!ionized, "CX/null selection cannot create ionization products");
            const bool exchanged = projectile.x() == 0;
            check::vector(va, exchanged ? stationary : incident, 0, "CX/null projectile velocity");
            check::vector(vb, exchanged ? incident : stationary, 0, "CX/null target velocity");
            if (replay == 0) {
                first_run.push_back(exchanged);
                reactions += exchanged;
            } else {
                check::equal(exchanged, first_run[trial], "seed must reproduce CX/null channel sequence");
            }
        }
    }
    const double sigma = std::sqrt(trials * 0.5 * 0.5);
    check::near((reactions - trials * 0.5) / sigma, 0, 6, "CX/null selection outside six-sigma interval");
}

// A fresh projectile meets a stationary neutral reservoir in every trial.
// The number of reactions must follow Binomial(N, 1-exp(-nu*dt)).
// This checks the absolute rate, probability sampler and channel selection.
void reaction_probability(bool electron, CollisionScheme scheme, ElectronIonizationModel energy_model) {
    constexpr int trials = 32768;
    constexpr double probability = 0.2;
    constexpr double density_cm3 = 1e13;
    CollisionProcessOptions options;
    options.electron_ionization = electron;
    options.proton_ionization = false;
    options.proton_charge_exchange = !electron;
    ColliderWithNeutrals collider(density_cm3, scheme, options);
    collider.rng.seed(1729);
    set_electron_ionization_model(energy_model);
    const double mass = electron ? 1 : 1836;
    const double energy = (electron ? 0.05 : 1.0) / 511.0;   // 50 eV or 1 keV, in m_e*c^2
    // Reference values of the implemented hydrogen fits at the specified energies,
    // in cm^2; these are regression points, not experimental cross-section data.
    const double cross_section_cm2 = electron ? 6.6255815613781293e-17 : 1.707817171554461e-15;
    const Vector3R incident(std::sqrt(2 * energy / mass), 0, 0), stationary(0.0, 0.0, 0.0);
    const double length_unit_cm = SGS::c / SGS::get_plasma_freq(density_cm3);
    const double expected_frequency = incident.norm() * cross_section_cm2 * density_cm3 * length_unit_cm;
    const double actual_frequency = collider.total_collision_frequency(incident, stationary, mass, 1);
    check::near(actual_frequency / expected_frequency, 1, 1e-12, "absolute hydrogen reaction frequency");
    const double dt = -std::log(1 - probability) / expected_frequency;
    int reacted = 0;
    double max_energy_error = 0;
    for (int trial = 0; trial < trials; ++trial) {
        auto projectile = incident, target = stationary;
        auto [ionized, secondary_electron, product_ion] = collider.collision_with_neutral(
            projectile, target, mass, electron ? 1837 : 1836, 1, dt, expected_frequency);
        if (electron) {
            if (!ionized) {
                check::vector(projectile, incident, 0, "no reaction must preserve projectile velocity");
                continue;
            }
            ++reacted;
            // The implemented electron-impact model neglects ion recoil:
            // two outgoing electrons share E_incident - I_H, the ion is at rest.
            const double outgoing_energy = 0.5 * (projectile.squaredNorm() + secondary_electron.squaredNorm());
            const double error = std::abs(outgoing_energy + THRESHOLD_IONIZATION - energy);
            check::near(error, 0, 1e-15, "ionization electron energy plus hydrogen binding energy");
            check::vector(product_ion, stationary, 0, "stationary ion in the no-recoil model");
            max_energy_error = std::max(max_energy_error, error);
        } else {
            check::require(!ionized, "charge exchange must not create an electron-ion pair");
            if ((projectile - incident).norm() == 0)
                continue;
            ++reacted;
            check::vector(projectile, stationary, 0, "CX projectile takes neutral velocity");
            check::vector(target, incident, 0, "CX neutral takes projectile velocity");
        }
    }
    const double expected_count = trials * probability;
    const double standard_deviation = std::sqrt(trials * probability * (1 - probability));
    const double z_score = (reacted - expected_count) / standard_deviation;
    std::cout << (electron ? "electron_ionization" : "charge_exchange")
              << (scheme == CollisionScheme::PHYSICAL_ONLY ? " physical_only" : " null_collision")
              << " model=" << (energy_model == ElectronIonizationModel::Opal ? "Opal" : "Approximation")
              << " N=" << trials << " reactions=" << reacted << " expected=" << expected_count << " z=" << z_score
              << " (limit 6 sigma), max energy error=" << max_energy_error << '\n';
    check::near(z_score, 0, 6, "reaction count outside the six-sigma binomial interval");
}

// Proton beam through neutrals at rest, charge exchange only. With beam weight w_b and neutral weight w_n, one
// exchange of a beam macroparticle moves w_b physical particles; w_b/w_n neutral macroparticles (integer part
// plus a Bernoulli remainder) take the beam velocity. Physical momentum is conserved on average, with the
// remainder variance p(1-p) (w_n v)^2 per exchange and the noise of hitting already exchanged neutrals.
void charge_exchange_with_unequal_weights() {
    for (const double ratio : {4.0, 0.25, 2.5}) {
        Domain domain;
        domain.init(Vector3I(1, 1, 1), Vector3R(1, 1, 1));
        const int neutral_count = 40000;
        const int beam_count = static_cast<int>(std::lround(neutral_count / ratio / 20));   // beam 20x less dense
        Species species;
        species.emplace("Beam",
                        std::make_unique<ParticlesArray>(check::species("Beam", 1, 1836, 0.05, beam_count), domain));
        species.emplace("H", std::make_unique<ParticlesArray>(check::species("H", 0, 1836, 1, neutral_count), domain));
        species.emplace("e", std::make_unique<ParticlesArray>(check::species("e", -1, 1, 1, 1), domain));
        species.emplace("p", std::make_unique<ParticlesArray>(check::species("p", 1, 1836, 1, 1), domain));
        const Vector3R beam_velocity(std::sqrt(2 * (1.0 / 511) / 1836), 0, 0);   // 1 keV protons
        for (int i = 0; i < beam_count; ++i)
            species.at("Beam")->add_particle(Particle(Vector3R(0.5, 0.5, 0.5), beam_velocity));
        for (int i = 0; i < neutral_count; ++i)
            species.at("H")->add_particle(Particle(Vector3R(0.5, 0.5, 0.5), Vector3R::Zero()));
        const double w_beam = species.at("Beam")->mpw(), w_neutral = species.at("H")->mpw();

        NeutralCollisionConfig config;
        config.charged_species = "Beam";
        config.neutral_species = "H";
        config.electron_product = "e";
        config.ion_product = "p";
        config.electron_ionization = false;
        config.proton_ionization = false;
        config.proton_charge_exchange = true;
        CollisionProcessOptions options;
        options.electron_ionization = false;
        options.proton_ionization = false;
        options.proton_charge_exchange = true;
        ColliderWithNeutrals rate(1e13, CollisionScheme::PHYSICAL_ONLY, options);
        // Exchange probability 0.3 per beam particle in one step.
        const double frequency = rate.total_collision_frequency(beam_velocity, Vector3R::Zero(), 1836, 1.0);
        NeutralCollisionOperator collisions(config, 1e13);
        collisions.apply(species, -std::log(0.7) / frequency);

        int exchanged = 0, neutrals_hit = 0;
        for (int cell = 0; cell < species.at("Beam")->size(); ++cell) {
            for (const Particle& p : (*species.at("Beam"))(cell)) exchanged += p.velocity.x() == 0;
            for (const Particle& p : (*species.at("H"))(cell)) neutrals_hit += p.velocity.x() != 0;
        }
        const double momentum_change = (neutrals_hit * w_neutral - exchanged * w_beam) * beam_velocity.x();
        // Sampling noise: the Bernoulli remainder of w_b/w_n per exchange, and hits on neutrals that already took
        // the beam velocity earlier in the step (fraction f), which change nothing but shift both counts.
        const double remainder = ratio - std::floor(ratio);
        const double fast = static_cast<double>(neutrals_hit) / neutral_count;
        const double sigma =
            std::sqrt(exchanged * (remainder * (1 - remainder) + ratio * fast)) * w_neutral * beam_velocity.x() +
            std::sqrt(exchanged * fast) * w_beam * beam_velocity.x();
        const double scale = w_beam * beam_count * beam_velocity.x();
        std::cout << "CX weights beam/neutral=" << ratio << ": beam exchanged " << exchanged << "/" << beam_count
                  << ", neutrals hit " << neutrals_hit << ", momentum change " << momentum_change / scale
                  << " of the beam momentum (sampling sigma " << sigma / scale << ")\n";
        check::require(std::abs(exchanged - 0.3 * beam_count) < 6 * std::sqrt(0.21 * beam_count),
                       "charge exchange count outside six sigma");
        check::require(std::abs(momentum_change) <= 5 * sigma + 1e-12 * scale,
                       "charge exchange with unequal weights does not conserve momentum on average");
    }
}

int main() {
    return check::run([] {
        std::cout << std::scientific << std::setprecision(8);
        zero_limits_and_charge_exchange();
        seeded_null_channel_selection();
        const auto original_model = get_electron_ionization_model();
        for (CollisionScheme scheme : {CollisionScheme::PHYSICAL_ONLY, CollisionScheme::NULL_COLLISION}) {
            reaction_probability(false, scheme, ElectronIonizationModel::Approximation);
            for (auto model : {ElectronIonizationModel::Approximation, ElectronIonizationModel::Opal})
                reaction_probability(true, scheme, model);
        }
        set_electron_ionization_model(original_model);
        charge_exchange_with_unequal_weights();
        std::cout << "PASS: neutral collisions (absolute rate, probability, CX, ionization energy)\n";
    });
}
