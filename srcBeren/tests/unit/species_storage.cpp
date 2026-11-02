#include <limits>
#include <map>
#include <type_traits>

#include "test_support.h"

static_assert(!std::is_copy_constructible_v<ParticleStorage>);
static_assert(!std::is_constructible_v<ParticleStorage, Domain&&>);
static_assert(!std::is_constructible_v<ParticleStorage, const Domain&&>);
static_assert(std::is_same_v<decltype(std::declval<ParticlesArray&>().parameters()), const SpeciesParameters&>);

namespace {
void moving_a_species_preserves_ownership_and_step_state() {
    Domain domain;
    domain.init(Vector3I(3, 3, 3), Vector3R(0.7, 1.1, 0.9));
    ParticlesArray species(check::species("beam", -2, 3, 2, 10), domain);
    const Vector3R position(0.8, 1.2, 1.0);
    species.add_particles({Particle(position, Vector3R(1, 2, 3)), Particle(position, Vector3R(4, 5, 6)),
                           Particle(position, Vector3R(7, 8, 9))});
    species.prepare();
    species.currentOnGrid.data().setConstant(0.75);
    species.diag.add_loss(Face::ZMIN, 0.125);
    const Vector3I index = domain.get_cell_index(position);
    const int cell = domain.sind(index.x(), index.y(), index.z());
    const Particle* addresses = species.cells()(cell).data();
    const double energy = species.kineticEnergy;

    ParticlesArray moved(std::move(species));
    check::require(&moved.get_domain() == &domain, "moving species keeps the borrowed domain");
    check::require(moved.cells()(cell).data() == addresses, "moving species transfers existing cell allocations");
    check::require(moved.name() == "beam", "moving species preserves its name");
    check::near(moved.mpw(), 0.2, 0, "moving species preserves weight");
    check::near(moved.kineticEnergy, energy, 0, "moving species preserves beginning-of-step energy");
    check::near(moved.diag.total_lost_energy(), 0.125, 0, "moving species preserves step events");
    check::near(moved.currentOnGrid[0], 0.75, 0, "moving species preserves grid current");
    moved.swap_and_pop_particle(cell, 1);
    check::equal(moved.get_total_num_of_particles(), 2, "swap-and-pop removes exactly one particle");
    check::vector(moved.cells()(cell)[1].velocity, Vector3R(7, 8, 9), 0, "last particle fills the removed slot");
    moved.swap_and_pop_particle(index.x(), index.y(), index.z(), 0);
    check::equal(moved.get_total_num_of_particles(), 1, "three-index removal addresses the same cell");
    check::vector(moved.cells()(cell)[0].velocity, Vector3R(7, 8, 9), 0, "removal preserves the survivor");
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
        for (auto& particle : species.cells()(cell))
            particle.coord = expected.at(static_cast<int>(particle.velocity.x())).coord;
    species.update_cells(domain);
    check::equal(species.get_total_num_of_particles(), id, "redistribution preserves particle count");
    std::map<int, int> occurrences;
    for (int cell = 0; cell < species.size(); ++cell)
        for (const auto& particle : species.cells()(cell)) {
            const int particle_id = static_cast<int>(particle.velocity.x());
            const auto& reference = expected.at(particle_id);
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
    const auto settled = species.cells();
    species.update_cells(domain);
    for (int cell = 0; cell < species.size(); ++cell) {
        check::equal(species.cells()(cell).size(), settled(cell).size(), "settled cell size");
        for (size_t i = 0; i < settled(cell).size(); ++i)
            check::vector(species.cells()(cell)[i].velocity, settled(cell)[i].velocity, 0, "settled order");
    }
}

void species_parameters_preserve_normalization_and_validation() {
    Domain domain;
    domain.init(Vector3I(3, 3, 3), Vector3R(0.7, 1.1, 0.9));
    auto config = check::species("beam", -2, 3, 2, 10);
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
            auto invalid = config;
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
        auto invalid = config;
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
        check::scenario("Species parameters: normalized weight, neutral identity and invalid inputs",
                        species_parameters_preserve_normalization_and_validation);
        check::scenario("Neighbor redistribution: unique particles, correct cells and conserved energy",
                        redistribution_preserves_particles_and_cell_membership);
        check::scenario("Moving a species transfers storage and step state; swap-and-pop preserves survivors",
                        moving_a_species_preserves_ownership_and_step_state);
    });
}
