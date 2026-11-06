#include <omp.h>

#include "coupling/density.h"
#include "test_support.h"

namespace {
// Density is in code units: each particle contributes q*density/npc, with no
// division by cell volume. Neutrals contribute density/npc instead.
constexpr double macro_weight = 0.25;
constexpr double tolerance = 2e-13;

// Independent cardinal B-splines, evaluated over every stored node. Neither
// production stencil construction nor the blocked/reference deposit is used.
template <ShapeType ST>
long double analytical_weight(long double distance) {
    const long double r = std::abs(distance);
    if constexpr (ST == ShapeType::Linear) {
        return r < 1 ? 1 - r : 0;
    } else {
        if (r <= 0.5L)
            return 0.75L - r * r;
        if (r < 1.5L)
            return 0.5L * (1.5L - r) * (1.5L - r);
        return 0;
    }
}

std::vector<Vector3R> sample_endpoints() {
    // Sparse particles exercise empty blocks and the final partial two-cell
    // blocks. Pairs straddle integer and half-integer boundaries on every axis.
    return {{0, 0, 0},
            {5, 6, 7},
            {0, 2.2, 3.3},
            {5, 2.2, 3.3},
            {1.2, 0, 3.3},
            {1.2, 6, 3.3},
            {1.2, 2.2, 0},
            {1.2, 2.2, 7},
            {1.999999, 2.999999, 3.999999},
            {2.000001, 3.000001, 4.000001},
            {2.499999, 3.499999, 4.499999},
            {2.500001, 3.500001, 4.500001},
            {4.83, 5.72, 6.61},
            {0.13, 0.27, 0.39}};
}

template <ShapeType ST>
void check_density_nodes(int ghost, int threads, bool neutral, const std::vector<Vector3R>& normalized_endpoints,
                         const std::string& sample) {
    Domain domain;
    domain.init(Vector3I(5, 6, 7), Vector3R(0.5, 1.25, 2), ghost);
    ParticlesArray species(check::species("density_test", neutral ? 0 : -2, 3, 2, 8), domain);
    std::vector<Vector3R> positions;
    for (const auto& u : normalized_endpoints) {
        const Vector3R position(u.x() * 0.5, u.y() * 1.25, u.z() * 2);
        positions.push_back(position);
        species.add_particle(Particle(position, Vector3R(0.0)));
    }
    check::near(species.mpw(), macro_weight, 0, "fixture macro weight");
    const double particle_weight = macro_weight * (neutral ? 1 : -2);
    const std::string context = sample + ", ghost=" + std::to_string(ghost) + ", threads=" + std::to_string(threads) +
                                (neutral ? ", neutral" : ", charged");
    omp_set_num_threads(threads);
    int actual_threads = 0;
#pragma omp parallel
    {
#pragma omp single
        actual_threads = omp_get_num_threads();
    }
    check::require(actual_threads == threads, context + " requested OpenMP team must be available");
    Field3d density(domain.size(), 1);
    for (bool reference : {false, true}) {
        // Each call must overwrite all nodes, including untouched ghosts.
        density.data().setConstant(17);
        if (reference)
            coupling::deposit_density_reference<ST>(species, density);
        else
            coupling::deposit_density<ST>(species, density);
        const std::string kernel = context + (reference ? ", reference" : ", blocked");
        long double sum = 0;
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k) {
                    long double expected = 0;
                    for (const auto& position : positions) {
                        // Recover normalized coordinates from the actual inserted
                        // physical point so rounding at a cell boundary is shared.
                        const long double ux = position.x() / domain.cell_size(0);
                        const long double uy = position.y() / domain.cell_size(1);
                        const long double uz = position.z() / domain.cell_size(2);
                        expected += particle_weight * analytical_weight<ST>(i - ghost - ux) *
                                    analytical_weight<ST>(j - ghost - uy) * analytical_weight<ST>(k - ghost - uz);
                    }
                    const double actual = density(i, j, k, 0);
                    check::near(actual, double(expected), tolerance,
                                kernel + " node " + std::to_string(i - ghost) + "," + std::to_string(j - ghost) + "," +
                                    std::to_string(k - ghost));
                    sum += actual;
                }
        check::near(double(sum), particle_weight * positions.size(), tolerance, kernel + " total charge/number");
    }
}

template <ShapeType ST>
void density_has_analytical_node_values() {
    const int previous_threads = omp_get_max_threads();
    const int previous_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    for (int ghost : {1, 2})
        for (int threads : {1, 4})
            for (bool neutral : {false, true}) {
                check_density_nodes<ST>(ghost, threads, neutral, {}, "empty species overwrites field");
                check_density_nodes<ST>(ghost, threads, neutral, {{1.25, 2.5, 3.75}}, "single particle polynomial");
                check_density_nodes<ST>(ghost, threads, neutral, sample_endpoints(), "wall and crossing endpoints");
                // A populated ghost bin must not disappear. Both shapes have
                // their full nonzero support inside even the one-node halo.
                check_density_nodes<ST>(ghost, threads, neutral, {{-0.2, 0.1, -0.1}}, "negative ghost-bin endpoint");
                if constexpr (ST == ShapeType::Linear) {
                    // At an exact outer halo node the Linear shape has only one
                    // nonzero weight. This populates the final partial bin block
                    // without putting nonzero support outside allocated storage.
                    check_density_nodes<ST>(ghost, threads, neutral, {{5.0 + ghost, 6.0 + ghost, 7.0 + ghost}},
                                            "populated partial block at outer halo node");
                }
            }
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("Linear density matches analytical node weights and conserves charge/number",
                        density_has_analytical_node_values<ShapeType::Linear>);
        check::scenario("Quadratic density matches analytical node weights and conserves charge/number",
                        density_has_analytical_node_values<ShapeType::Quadratic>);
    });
}
