#include <algorithm>
#include <limits>
#include <map>
#include <type_traits>

#include "particles/particle_sources.h"
#include "test_support.h"
#include "utils/util.h"

static_assert(std::is_same_v<decltype(std::declval<ParticlesArray&>().parameters()), const SpeciesParameters&>);

namespace {
void generation_preserves_sampling_energy_and_rng(int count, bool bound, bool neutral) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    ParticlesArray species(check::species("beam", neutral ? 0 : -1, 2, 2, 4), domain);
    const Vector3R center(3.8, 2, 2), extent(0.3, 0.5, 0.5), velocity(0.5, 0.25, 0.125);
    const double dt = 0.25;
    InitialDistribution distribution(std::make_shared<RectangleDistribution>(center, extent),
                                     std::make_shared<GaussianVelocity>(velocity, Vector3R(0, 0, 0)),
                                     bound ? "injection_bound" : "initial", count, 2, 0.5);
    LehmerEngine space(123), momentum(456);
    const InjectedParticles injected =
        add_particles_from_distribution(species, distribution, space, momentum, domain, dt, true);
    std::vector<Vector3R> expected, actual;
    // Independently reconstruct the uniform coordinates from the specified
    // per-particle stream, then apply the boundary-injection drift and clipping.
    for (int i = 0; i < count; ++i) {
        LehmerEngine stream(123);
        stream.discard(32 * i);
        std::uniform_real_distribution<double> uniform(0, 1);
        const double x = center.x() + extent.x() * (1 - 2 * uniform(stream));
        const double y = center.y() + extent.y() * (1 - 2 * uniform(stream));
        const double z = center.z() + extent.z() * (1 - 2 * uniform(stream));
        const Vector3R position = Vector3R(x, y, z) + (bound ? dt * velocity : Vector3R(0, 0, 0));
        if (position.x() >= 0 && position.x() <= 4 && position.y() >= 0 && position.y() <= 4 && position.z() >= 0 &&
            position.z() <= 4)
            expected.push_back(position);
    }
    for (int cell = 0; cell < species.size(); ++cell)
        for (const Particle& particle : species(cell)) {
            actual.push_back(particle.coord);
            check::vector(particle.velocity, velocity, 0, "prescribed injection velocity");
        }
    const auto less = [](const Vector3R& a, const Vector3R& b) {
        if (a.x() != b.x())
            return a.x() < b.x();
        if (a.y() != b.y())
            return a.y() < b.y();
        return a.z() < b.z();
    };
    std::sort(actual.begin(), actual.end(), less);
    std::sort(expected.begin(), expected.end(), less);
    check::equal(actual.size(), expected.size(), "clipping preserves exactly the accepted samples");
    check::require(!actual.empty() && actual.size() < static_cast<size_t>(count), "exercise acceptance and rejection");
    for (size_t i = 0; i < actual.size(); ++i)
        check::vector(actual[i], expected[i], 0, "coordinate stream and bound-injection drift");
    check::equal(injected.count, static_cast<int>(actual.size()), "generation result counts accepted particles");
    check::near(injected.energy, actual.size() * 0.5 * 2 * 0.5 * velocity.squaredNorm(), 1e-12,
                "injection energy includes only accepted particles");
    LehmerEngine expected_space(123), expected_momentum(456);
    expected_space.discard(32 * count);
    expected_momentum.discard(32 * count);
    check::require(space.state == expected_space.state && momentum.state == expected_momentum.state,
                   "both RNG streams advance for rejected as well as accepted particles");
}

void fractional_injection_matches_rate_and_accepted_particles() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    for (double rate : {0.25, 0.75, 1.25}) {
        for (bool clipped : {false, true}) {
            ParticlesArray species(check::species("beam", -1, 2, 2, 4), domain);
            const Vector3R center(clipped ? 3.8 : 2, 2, 2), extent(clipped ? 0.3 : 0.5, 0.5, 0.5);
            std::vector<std::unique_ptr<IDistribution>> distributions;
            distributions.push_back(std::make_unique<InjectionDistribution>(
                std::make_shared<RectangleDistribution>(center, extent),
                std::make_shared<GaussianVelocity>(Vector3R(0.5, 0, 0), Vector3R(0, 0, 0)),
                clipped ? "injection_bound" : "injection", rate, 2, 0.5));
            int total_accepted = 0;
            for (int step = 1; step <= 12; ++step) {
                species.begin_step();
                const double energy = inject_particles_step(species, distributions, step, domain, 0.25, false);
                const int requested = static_cast<int>(std::floor(step * rate) - std::floor((step - 1) * rate));
                int accepted = 0;
                for (int i = 0; i < requested; ++i) {
                    LehmerEngine stream(13 * 3 * step + hash("beam", 1000003));
                    stream.discard(32 * i);
                    std::uniform_real_distribution<double> uniform(0, 1);
                    const double x = center.x() + extent.x() * (1 - 2 * uniform(stream)) + (clipped ? 0.125 : 0);
                    accepted += x >= 0 && x <= 4;
                }
                total_accepted += accepted;
                check::equal(species.step_state.injection_count, accepted, "step count includes only accepted samples");
                check::equal(species.get_total_num_of_particles(), total_accepted, "fractional injection population");
                check::near(energy, accepted * 0.125, 0, "accepted injection energy");
                if (!clipped)
                    check::equal(total_accepted, static_cast<int>(std::floor(step * rate)),
                                 "integrated fractional rate");
            }
        }
    }
}

void moving_a_species_preserves_ownership_and_step_state() {
    Domain domain;
    domain.init(Vector3I(3, 3, 3), Vector3R(0.7, 1.1, 0.9));
    ParticlesArray species(check::species("beam", -2, 3, 2, 10), domain);
    const Vector3R position(0.8, 1.2, 1.0);
    species.add_particles({Particle(position, Vector3R(1, 2, 3)), Particle(position, Vector3R(4, 5, 6)),
                           Particle(position, Vector3R(7, 8, 9))});
    species.begin_step();
    species.currentOnGrid.data().setConstant(0.75);
    species.step_state.add_loss(Face::ZMIN, 0.125);
    const Vector3I index = domain.get_cell_index(position);
    const int cell = domain.sind(index.x(), index.y(), index.z());
    const Particle* addresses = species(cell).data();
    const double energy = species.step_state.kinetic_energy_before_step;

    ParticlesArray moved(std::move(species));
    check::require(&moved.get_domain() == &domain, "moving species keeps the borrowed domain");
    check::require(moved(cell).data() == addresses, "moving species transfers existing cell allocations");
    check::require(moved.name() == "beam", "moving species preserves its name");
    check::near(moved.mpw(), 0.2, 0, "moving species preserves weight");
    check::near(moved.step_state.kinetic_energy_before_step, energy, 0,
                "moving species preserves beginning-of-step energy");
    check::near(moved.step_state.total_lost_energy(), 0.125, 0, "moving species preserves step events");
    check::near(moved.currentOnGrid[0], 0.75, 0, "moving species preserves grid current");
    moved.swap_and_pop_particle(cell, 1);
    check::equal(moved.get_total_num_of_particles(), 2, "swap-and-pop removes exactly one particle");
    check::vector(moved(cell)[1].velocity, Vector3R(7, 8, 9), 0, "last particle fills the removed slot");
    moved.swap_and_pop_particle(index.x(), index.y(), index.z(), 0);
    check::equal(moved.get_total_num_of_particles(), 1, "three-index removal addresses the same cell");
    check::vector(moved(cell)[0].velocity, Vector3R(7, 8, 9), 0, "removal preserves the survivor");
}

void redistribution_preserves_particles_and_cell_membership() {
    Domain domain;
    domain.init(Vector3I(8, 8, 8), Vector3R(0.7, 1.1, 0.9));
    ParticlesArray species(check::species("beam", -2, 3, 2, 10), domain);
    std::map<int, Particle> expected;
    int id = 0;
    // All 27 offsets, including stationary particles; both directions of each
    // axis and all coloring classes are populated. Every displacement is < dx.
    for (int x = 2; x < 6; ++x)
        for (int y = 2; y < 6; ++y)
            for (int z = 2; z < 6; ++z)
                for (int dx = -1; dx <= 1; ++dx)
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dz = -1; dz <= 1; ++dz) {
                            const Vector3R position((x + 0.5) * 0.7, (y + 0.5) * 1.1, (z + 0.5) * 0.9);
                            // Velocity identifies the particle without modifying Particle's layout.
                            ++id;
                            const Vector3R velocity(id, -0.5 * id, 0.25 * id);
                            species.add_particle(Particle(position, velocity));
                            expected.emplace(
                                id, Particle(position + Vector3R(dx * 0.6 * 0.7, dy * 0.6 * 1.1, dz * 0.6 * 0.9),
                                             velocity));
                        }
    const double energy_before = species.get_kinetic_energy();
    for (int cell = 0; cell < species.size(); ++cell)
        for (Particle& particle : species(cell))
            particle.coord = expected.at(static_cast<int>(particle.velocity.x())).coord;
    species.update_cells(domain);
    check::equal(species.get_total_num_of_particles(), id, "redistribution preserves particle count");
    std::map<int, int> occurrences;
    for (int cell = 0; cell < species.size(); ++cell)
        for (const Particle& particle : species(cell)) {
            const int particle_id = static_cast<int>(particle.velocity.x());
            const Particle& reference = expected.at(particle_id);
            ++occurrences[particle_id];
            check::vector(particle.coord, reference.coord, 0, "redistribution changes no coordinate");
            check::vector(particle.velocity, reference.velocity, 0, "redistribution changes no velocity");
            const Vector3I index = domain.get_cell_index(reference.coord);
            check::equal(cell, domain.sind(index.x(), index.y(), index.z()),
                         "particle belongs to its destination cell");
        }
    for (const auto& [particle_id, particle] : expected)
        check::equal(occurrences[particle_id], 1, "each particle occurs exactly once");
    check::near(species.get_kinetic_energy(), energy_before, 1e-6, "redistribution conserves kinetic energy");
    // A second pass must leave the order within every cell unchanged.
    const Array3D<std::vector<Particle>> settled = species.cells();
    species.update_cells(domain);
    for (int cell = 0; cell < species.size(); ++cell) {
        check::equal(species(cell).size(), settled(cell).size(), "settled cell size");
        for (size_t i = 0; i < settled(cell).size(); ++i)
            check::vector(species(cell)[i].velocity, settled(cell)[i].velocity, 0, "settled order");
    }
}

void species_parameters_preserve_normalization_and_validation() {
    Domain domain;
    domain.init(Vector3I(3, 3, 3), Vector3R(0.7, 1.1, 0.9));
    nlohmann::json config = check::species("beam", -2, 3, 2, 10);
    ParticlesArray beam(config, domain);
    check::near(beam.charge(), -2, 0, "configured charge");
    check::near(beam.mass(), 3, 0, "configured mass");
    // mpw is density / particles per cell, without a cell-volume factor.
    check::near(beam.mpw(), 0.2, 0, "macroparticle weight normalization");
    check::require(beam.name() == "beam" && !beam.is_neutral(), "charged species identity");
    ParticlesArray neutral(check::species("gas", 0), domain);
    ParticlesArray named_neutral(check::species("Neutrals", 1), domain);
    check::require(neutral.is_neutral() && named_neutral.is_neutral(), "neutral classification contract");
    ParticlesArray empty_density(check::species("empty", -1, 1, 0), domain);
    check::near(empty_density.mpw(), 0, 0, "zero density remains supported");
    for (const std::string key : {"Charge", "Mass", "Density"})
        for (double value : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
            nlohmann::json invalid = config;
            invalid[key] = value;
            bool rejected = false;
            try {
                ParticlesArray bad(invalid, domain);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            check::require(rejected, "nonfinite " + key + " must be rejected");
        }
    for (const auto& [key, value] :
         std::map<std::string, double>{{"Mass", 0}, {"Density", -1}, {"NumPartPerCell", 0}}) {
        nlohmann::json invalid = config;
        invalid[key] = value;
        bool rejected = false;
        try {
            ParticlesArray bad(invalid, domain);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check::require(rejected, "invalid " + key + " must be rejected");
    }
}
}   // namespace

int main() {
    return check::run([] {
        for (int count : {31, 15017})
            for (bool bound : {false, true})
                for (bool neutral : {false, true})
                    check::scenario(
                        "Generation: RNG streams, clipping and accepted energy; count=" + std::to_string(count) +
                            (bound ? " bound" : " initial") + (neutral ? " neutral" : " charged"),
                        [&] { generation_preserves_sampling_energy_and_rng(count, bound, neutral); });
        check::scenario("Fractional injection: one rate sample per step, accepted count and energy",
                        fractional_injection_matches_rate_and_accepted_particles);
        check::scenario("Species parameters: normalized weight, neutral identity and invalid inputs",
                        species_parameters_preserve_normalization_and_validation);
        check::scenario("Neighbor redistribution: unique particles, correct cells and conserved energy",
                        redistribution_preserves_particles_and_cell_membership);
        check::scenario("Moving a species transfers storage and step state; swap-and-pop preserves survivors",
                        moving_a_species_preserves_ownership_and_step_state);
    });
}
