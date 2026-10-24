#include "collisions_with_neutrals.h"
#include "test_support.h"

void test_neutrals() {
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
int main() {
    return check::run(test_neutrals);
}
