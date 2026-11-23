#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "particles/particle_sources.h"
#include "test_support.h"

namespace {
void van_der_corput_values() {
    check::near(van_der_corput(1, 2), 0.5, 0.0, "base-2 first value");
    check::near(van_der_corput(2, 2), 0.25, 0.0, "base-2 second value");
    check::near(van_der_corput(3, 2), 0.75, 0.0, "base-2 third value");
    check::near(van_der_corput(1, 3), 1.0 / 3.0, 1e-15, "base-3 first value");
}

void halton_positions_are_deterministic_and_distinct() {
    HaltonDistribution distribution(Vector3R(0, 0, 0), Vector3R(2, 2, 2), Vector3R(0, 0, 0),
                                                  Vector3R(0, 0, 0), 16, 1.0, 1.0);
    check::require(distribution.get_count_to_inject() == 16, "Halton distribution count");
    LehmerEngine first(1), second(1);
    std::vector<Vector3R> positions;
    for (std::size_t i = 0; i < 16; ++i) {
        Vector3R p1, v1, p2, v2;
        distribution.sample_particle(i, first, first, p1, v1);
        distribution.sample_particle(i, second, second, p2, v2);
        check::vector(p1, p2, 0.0, "Halton position is deterministic");
        check::vector(v1, Vector3R(0, 0, 0), 0.0, "zero sigma gives zero velocity");
        for (int d = 0; d < 3; ++d)
            check::require(std::abs(p1[d]) <= 2.0, "Halton position stays inside the rectangle");
        for (std::size_t j = 0; j < positions.size(); ++j)
            check::require((p1 - positions[j]).length() > 1e-14, "Halton positions are distinct");
        positions.push_back(p1);
    }
}

void quiet_start_pairs() {
    QuietStartDistribution distribution(Vector3R(0, 0, 0), Vector3R(1, 1, 1), Vector3R(0, 0, 0),
                                        Vector3R(0.1, 0.2, 0.3), 8, 2, 1.0, 1.0);
    check::require(distribution.get_count_to_inject() == 8, "quiet-start count");
    LehmerEngine engine(3);
    for (std::size_t pair = 0; pair < 4; ++pair) {
        Vector3R positive_position, positive_velocity, negative_position, negative_velocity;
        distribution.sample_particle(2 * pair, engine, engine, positive_position, positive_velocity);
        distribution.sample_particle(2 * pair + 1, engine, engine, negative_position, negative_velocity);
        check::vector(positive_position, negative_position, 0.0, "quiet pair shares a coordinate");
        for (int d = 0; d < 3; ++d)
            check::near(positive_velocity[d], -negative_velocity[d], 1e-13, "quiet pair has antithetic velocity");
        check::near(positive_position.z(), -0.75 + 0.5 * pair, 1e-13, "quiet stratification along z");
    }
}

// n(s) ~ 1 + a sin(2 pi m s + phi) along the stratified axis: the markers
// follow the modulated cumulative fraction, and <sin(2 pi m s + phi)> = a / 2.
void quiet_start_density_perturbation() {
    const SinePerturbation perturbation{0.2, 2, 0.3, 0.05};
    const int count = 4000;
    QuietStartDistribution distribution(Vector3R(0, 0, 4), Vector3R(1, 1, 4), Vector3R(0, 0, 0), Vector3R(0, 0, 0),
                                        count, 2, 1.0, 1.0, perturbation);
    LehmerEngine engine(5);
    double mean_sine = 0, worst_cdf = 0;
    const double k = 2 * M_PI * perturbation.mode;
    for (int i = 0; i < count; i += 2) {
        Vector3R position, velocity;
        distribution.sample_particle(i, engine, engine, position, velocity);
        const double s = position.z() / 8.0;   // fraction across the box z in [0, 8)
        // sigma = 0: the velocity is the local mean flow.
        check::near(velocity.z(), 0.05 * std::sin(k * s + perturbation.phase), 1e-15, "mean flow of the wave");
        mean_sine += 2.0 * std::sin(k * s + perturbation.phase) / count;
        const double cdf =
            s + perturbation.amplitude / k * (std::cos(perturbation.phase) - std::cos(k * s + perturbation.phase));
        worst_cdf = std::max(worst_cdf, std::abs(cdf - (i / 2 + 0.5) / (count / 2)));
    }
    check::near(worst_cdf, 0, 1e-13, "markers sit at the modulated cumulative fraction");
    check::near(mean_sine, 0.5 * perturbation.amplitude, 1e-5, "first harmonic of the loaded density");
}

void halton_velocity_is_finite() {
    HaltonDistribution distribution(Vector3R(0, 0, 0), Vector3R(1, 1, 1), Vector3R(1, 2, 3),
                                                  Vector3R(0.1, 0.2, 0.3), 32, 1.0, 1.0);
    LehmerEngine engine(7);
    for (std::size_t i = 0; i < 32; ++i) {
        Vector3R position, velocity;
        distribution.sample_particle(i, engine, engine, position, velocity);
        bool finite = true;
        for (int d = 0; d < 3; ++d)
            finite = finite && std::isfinite(position[d]) && std::isfinite(velocity[d]);
        check::require(finite, "Halton sample is finite");
    }
}
}  // namespace

int main() {
    return check::run([] {
        check::scenario("Van der Corput sequence", van_der_corput_values);
        check::scenario("Halton deterministic positions", halton_positions_are_deterministic_and_distinct);
        check::scenario("Halton finite velocities", halton_velocity_is_finite);
        check::scenario("Halton quiet-start pairs", quiet_start_pairs);
        check::scenario("Halton quiet-start density perturbation", quiet_start_density_perturbation);
    });
}
