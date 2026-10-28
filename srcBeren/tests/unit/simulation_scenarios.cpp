#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>

#include "simulation_ecsim_corr.h"
#include "test_support.h"

namespace {
constexpr int cells = 5;
constexpr int steps = 20;
constexpr double dt = 0.2;
constexpr double energy_tolerance = 1e-7;

nlohmann::json settings(const std::string& scheme) {
    return {{"Scheme", scheme},
            {"Dt", dt},
            {"NumCellsX", cells},
            {"NumCellsY", cells},
            {"NumCellsZ", cells},
            {"Dx", 1},
            {"Dy", 1},
            {"Dz", 1},
            {"StartFromTime", 0},
            {"StartTimeStep", 0},
            {"LastTimestep", steps},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collider", "None"},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

// Observe every completed production step without reproducing its orchestration
// in the test. The observer reads fields/particles independently of energyConserve.
template <typename Scheme>
class PeriodicPlasma : public Scheme {
   public:
    using Scheme::Scheme;
    int samples = 0;
    double initial_energy = 0;
    double max_energy_error = 0;
    double max_field_energy = 0;
    double max_continuity_error = 0;
    Field3d previous_charge;
    std::ofstream history{"invariants.csv"};

    void make_diagnostic(int step) override {
        Scheme::make_diagnostic(step);
        const auto& domain = this->domain;
        Field3d charge(domain.size(), 1);
        this->collect_charge_density(charge);
        double kinetic = 0, field_energy = 0;
        // Integrate a single periodic copy: ghost nodes must not be counted twice.
        const int first = domain.grid.ghost_cells();
        for (int i = first; i < first + cells; ++i)
            for (int j = first; j < first + cells; ++j)
                for (int k = first; k < first + cells; ++k)
                    for (int d = 0; d < 3; ++d) {
                        const double e = this->fieldEn(i, j, k, d), b = this->fieldBn(i, j, k, d);
                        check::require(std::isfinite(e) && std::isfinite(b), "nonfinite periodic field");
                        field_energy += 0.5 * (e * e + b * b);   // cell volume = 1
                    }
        for (const auto& [name, species] : this->species) {
            check::near(species->get_total_num_of_particles(), 2 * cells * cells * cells, 0,
                        "periodicity must conserve each species particle count");
            double deposited_charge = 0;
            for (int i = first; i < first + cells; ++i)
                for (int j = first; j < first + cells; ++j)
                    for (int k = first; k < first + cells; ++k) deposited_charge += species->densityOnGrid(i, j, k, 0);
            // Each species carries nonzero charge; checking only their zero sum
            // would also accept a broken deposit returning zero everywhere.
            check::near(deposited_charge, species->charge * cells * cells * cells, 1e-10,
                        "deposited species charge differs from particle charge");
            // The diagnostic reads the stored density; it must equal a fresh
            // deposit of the current particles (a frozen/stale grid would still
            // integrate to the right total charge above). deposit_density()
            // applies the scheme shape and the density boundary condition.
            Field3d stored_density = species->densityOnGrid;
            this->deposit_density(*species);
            check::near((species->densityOnGrid - stored_density).norm(), 0, 1e-12,
                        "stored density must match a fresh deposit");
            species->densityOnGrid = stored_density;
            for (int cell = 0; cell < species->size(); ++cell)
                for (const auto& p : species->particlesData(cell)) {
                    for (int d = 0; d < 3; ++d) {
                        check::require(std::isfinite(p.velocity[d]) && std::isfinite(p.coord[d]) && p.coord[d] >= 0 &&
                                           p.coord[d] < cells,
                                       "particle must remain finite and inside the periodic box");
                    }
                    kinetic += 0.5 * species->mass() * species->mpw() * p.velocity.squaredNorm();
                }
        }
        const double total = kinetic + field_energy;
        if (step == 0) {
            previous_charge.resize(domain.size(), 1);
            initial_energy = total;
            check::require(initial_energy > 0, "periodic plasma needs nonzero thermal energy");
            history << "step,kinetic,field,total,relative_energy_error,max_continuity_error\n";
        }
        const double energy_error = std::abs(total - initial_energy) / initial_energy;
        check::near(energy_error, 0, energy_tolerance, "periodic total energy changed");
        if constexpr (std::is_same_v<Scheme, SimulationEcsimCorr>) {
            if (step > 0) {
                const Eigen::VectorXd div_current = this->mesh.divE * this->fieldJe.data();
                for (int i = first; i < first + cells; ++i)
                    for (int j = first; j < first + cells; ++j)
                        for (int k = first; k < first + cells; ++k) {
                            const double residual = (charge(i, j, k, 0) - previous_charge(i, j, k, 0)) / dt +
                                                    div_current[domain.sind(i, j, k)];
                            check::near(residual, 0, 1e-11, "periodic Corr local charge continuity");
                            max_continuity_error = std::max(max_continuity_error, std::abs(residual));
                        }
            }
        }
        previous_charge = charge;
        max_energy_error = std::max(max_energy_error, energy_error);
        max_field_energy = std::max(max_field_energy, field_energy);
        history << std::setprecision(17) << step << ',' << kinetic << ',' << field_energy << ',' << total << ','
                << energy_error << ',' << max_continuity_error << '\n';
        check::require(bool(history), "cannot write periodic invariant history");
        ++samples;
    }
};

template <typename Scheme>
void exercise(const std::string& scheme) {
    auto particles = nlohmann::json{{"particles", nlohmann::json::array({check::species("Electrons", -1, 1, 1, 2),
                                                                         check::species("Ions", 1, 1836, 1, 2)})}};
    char program[] = "simulation_scenarios";
    char* argv[] = {program};
    PeriodicPlasma<Scheme> simulation(settings(scheme), particles, 1, argv);
    simulation.init();
    // Antithetic thermal pairs suppress initial bulk current without making the
    // evolving deposit identically zero. Fixed seed, equal electron/ion T.
    std::mt19937 engine(1729);
    std::normal_distribution<double> thermal(0, 0.04);
    for (int i = 0; i < cells; ++i)
        for (int j = 0; j < cells; ++j)
            for (int k = 0; k < cells; ++k) {
                const Vector3R position(i + 0.98, j + 0.35, k + 0.65);
                const Vector3R velocity(thermal(engine), thermal(engine), thermal(engine));
                for (const auto& [name, species] : simulation.species)
                    for (int sign : {-1, 1})
                        species->add_particle(Particle(position, velocity * (sign / std::sqrt(species->mass()))));
            }
    for (auto& [name, species] : simulation.species) {
        simulation.deposit_density(*species);
        species->kineticEnergy = species->get_kinetic_energy();
    }
    simulation.calculate();
    simulation.finalize();
    check::near(simulation.samples, steps + 1, 0, "all production steps must be observed");
    check::require(simulation.max_field_energy > 1e-10 * simulation.initial_energy,
                   "periodic plasma must exercise a nonzero self-consistent field");
    std::cout << scheme << ": max relative energy error=" << simulation.max_energy_error
              << ", max field/initial energy=" << simulation.max_field_energy / simulation.initial_energy
              << ", max continuity error=" << simulation.max_continuity_error << '\n';
}
}   // namespace

int main() {
    return check::run([] {
        const auto root = std::filesystem::current_path();
        for (const std::string scheme : {"ecsim", "ecsim_corr"}) {
            const auto directory = root / scheme;
            std::filesystem::create_directories(directory);
            std::filesystem::current_path(directory);
            if (scheme == "ecsim")
                exercise<SimulationEcsim>(scheme);
            else
                exercise<SimulationEcsimCorr>(scheme);
            std::filesystem::current_path(root);
        }
    });
}
