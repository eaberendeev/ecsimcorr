#include <iomanip>

#include "collisions_with_neutrals.h"
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
    const Vector3R incident(std::sqrt(2 * energy / mass), 0, 0), stationary(0.0);
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

int main() {
    return check::run([] {
        std::cout << std::scientific << std::setprecision(8);
        zero_limits_and_charge_exchange();
        const auto original_model = get_electron_ionization_model();
        for (CollisionScheme scheme : {CollisionScheme::PHYSICAL_ONLY, CollisionScheme::NULL_COLLISION}) {
            reaction_probability(false, scheme, ElectronIonizationModel::Approximation);
            for (auto model : {ElectronIonizationModel::Approximation, ElectronIonizationModel::Opal})
                reaction_probability(true, scheme, model);
        }
        set_electron_ionization_model(original_model);
        std::cout << "PASS: neutral collisions (absolute rate, probability, CX, ionization energy)\n";
    });
}
