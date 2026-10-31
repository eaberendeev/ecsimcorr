#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>

#include "DensityMeasurements.h"
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
    std::vector<double> observed_state;
    std::vector<std::array<double, 3>> expected_continuity;
    double initial_energy = 0;
    double max_energy_error = 0;
    double max_field_energy = 0;
    double max_continuity_error = 0;
    Field3d previous_charge;
    std::ofstream history{"invariants.csv"};

    void make_diagnostic(int step) override {
        Scheme::make_diagnostic(step);
        for (const auto* field : {&this->fieldEn, &this->fieldBn})
            for (int i = 0; i < field->capacity(); ++i) observed_state.push_back((*field)[i]);
        for (const auto& [name, species] : this->species) {
            for (int i = 0; i < species->currentOnGrid.capacity(); ++i)
                observed_state.push_back(species->currentOnGrid[i]);
            for (int cell = 0; cell < species->size(); ++cell)
                for (const auto& particle : species->particlesData(cell))
                    for (int d = 0; d < 3; ++d) {
                        observed_state.push_back(particle.coord[d]);
                        observed_state.push_back(particle.velocity[d]);
                    }
        }
        const auto& domain = this->domain;
        const auto charge =
            measure_charge_density(this->species, domain, this->bc_handler, this->diagnostic_density_shape());
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
            const auto density = measure_density(*species, domain, this->bc_handler, this->diagnostic_density_shape());
            check::near(species->get_total_num_of_particles(), 2 * cells * cells * cells, 0,
                        "periodicity must conserve each species particle count");
            double deposited_charge = 0;
            for (int i = first; i < first + cells; ++i)
                for (int j = first; j < first + cells; ++j)
                    for (int k = first; k < first + cells; ++k) deposited_charge += density(i, j, k, 0);
            // Each species carries nonzero charge; checking only their zero sum
            // would also accept a broken deposit returning zero everywhere.
            check::near(deposited_charge, species->charge * cells * cells * cells, 1e-10,
                        "deposited species charge differs from particle charge");
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
                Eigen::VectorXd div_current(this->domain.total_size());
                spmv(this->maxwell.divE, this->fieldJe.data(), div_current);
                if (step % 3 == 0) {
                    const Eigen::VectorXd residual = (charge.data() - previous_charge.data()) / dt + div_current;
                    expected_continuity.push_back({step * dt, residual.norm(), residual.cwiseAbs().maxCoeff()});
                }
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
std::vector<double> exercise(const std::string& scheme, bool with_outputs) {
    auto particles = nlohmann::json{{"particles", nlohmann::json::array({check::species("Electrons", -1, 1, 1, 2),
                                                                         check::species("Ions", 1, 1836, 1, 2)})}};
    char program[] = "simulation_scenarios";
    char* argv[] = {program};
    auto config = settings(scheme);
    if (with_outputs)
        config["diagnostics"]["outputs"] = nlohmann::json::parse(R"([
            {"type":"energy_balance","interval":1},
            {"type":"pressure_2d","species":["Electrons"],"planes":{"z":[3]},"interval":2,"average_steps":3},
            {"type":"energy_spectrum","species":["Electrons"],"interval":2},
            {"type":"total_current_2d","planes":{"z":[3]},"interval":2},
            {"type":"probes","points":[[2,2,2]],"species":["Electrons"],"interval":2}
        ])");
    if constexpr (std::is_same_v<Scheme, SimulationEcsimCorr>)
        if (with_outputs)
            config["diagnostics"]["outputs"].push_back({{"type", "charge_continuity"}, {"interval", 3}});
    PeriodicPlasma<Scheme> simulation(config, particles, 1, argv);
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
        species->kineticEnergy = species->get_kinetic_energy();
    }
    simulation.calculate();
    simulation.finalize();
    if constexpr (std::is_same_v<Scheme, SimulationEcsimCorr>) {
        if (with_outputs) {
            std::ifstream file("continuity.txt");
            std::string header;
            std::getline(file, header);
            check::require(header == "Time L2 Linf", "missing continuity header");
            for (const auto& expected : simulation.expected_continuity)
                for (double value : expected) {
                    double actual;
                    check::require(bool(file >> actual), "missing selected continuity sample");
                    check::near(actual, value, 1e-12,
                                "scheduled residual differs from adjacent-step densities measured independently");
                }
            file >> std::ws;
            check::require(file.peek() == std::char_traits<char>::eof(), "unexpected continuity rows");
        }
    }
    check::near(simulation.samples, steps + 1, 0, "all production steps must be observed");
    check::require(simulation.max_field_energy > 1e-10 * simulation.initial_energy,
                   "periodic plasma must exercise a nonzero self-consistent field");
    std::cout << scheme << ": max relative energy error=" << simulation.max_energy_error
              << ", max field/initial energy=" << simulation.max_field_energy / simulation.initial_energy
              << ", max continuity error=" << simulation.max_continuity_error << '\n';
    return std::move(simulation.observed_state);
}
}   // namespace

int main() {
    return check::run([] {
        const auto root = std::filesystem::current_path();
        for (const std::string scheme : {"ecsim", "ecsim_corr"}) {
            std::vector<double> without_outputs;
            for (bool enabled : {false, true}) {
                const auto directory = root / (scheme + (enabled ? "_outputs" : ""));
                std::filesystem::create_directories(directory);
                std::filesystem::current_path(directory);
                const auto state = scheme == "ecsim" ? exercise<SimulationEcsim>(scheme, enabled)
                                                     : exercise<SimulationEcsimCorr>(scheme, enabled);
                if (enabled) {
                    check::require(state.size() == without_outputs.size(), "output changed observed state size");
                    for (size_t i = 0; i < state.size(); ++i)
                        check::near(state[i], without_outputs[i], 1e-12,
                                    "outputs changed a field, deposited current or particle during evolution");
                } else {
                    without_outputs = state;
                }
                std::filesystem::current_path(root);
            }
        }
    });
}
