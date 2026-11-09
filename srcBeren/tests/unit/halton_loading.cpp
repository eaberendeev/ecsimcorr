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
    });
}
