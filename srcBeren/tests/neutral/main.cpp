// Neutral-gas ionization by a beam, with equal and unequal macro weights.
//
// A beam (electrons at 1 keV or protons at 50 keV) ten times denser than the neutral gas ionizes it in one cell.
// The neutral density must follow d ln n_n / dt = -nu(t), nu(t) = sum over beam macroparticles of w_p sigma(v_p) v_p,
// evaluated from the actual beam state at every step, so a slowly changing beam does not bias the check. The
// cross sections themselves are pinned by neutral_collision_invariants; this experiment checks reaction sampling,
// the weight bookkeeping (neutrals removed, electrons and ions created) and the neutral density normalization.
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <random>
#include <string>
#include <vector>

#include "collisions/neutral/collisions_with_neutrals.h"
#include "collisions/neutral_collision.h"
#include "config/domain_config.h"
#include "grid/domain.h"
#include "particles/particle_array.h"

namespace {
constexpr double n0 = 1e13;
constexpr double ion_mass = 1836.0;

struct SpeciesWeight {
    double density;   // physical density in units of n0
    int count;        // macroparticles; weight = density / count (one cell of unit volume)
};

struct IonizationCase {
    std::string name;
    bool electron_beam;
    SpeciesWeight beam, neutrals, electrons, ions;   // electrons/ions: weight only (count = ppc reference)
};

nlohmann::json species_config(const std::string& name, double charge, double mass, double density, int ppc) {
    return {{"Name", name},       {"Charge", charge},      {"Mass", mass},
            {"Density", density}, {"NumPartPerCell", ppc}, {"distribution", nlohmann::json::array()}};
}

double physical(const ParticlesArray& sp) {
    double total = 0;
    for (int cell = 0; cell < sp.size(); ++cell) total += sp(cell).size();
    return total * sp.mpw();
}

int macroparticles(const ParticlesArray& sp) {
    int total = 0;
    for (int cell = 0; cell < sp.size(); ++cell) total += sp(cell).size();
    return total;
}

// Ionization rate of one neutral: sum over beam macroparticles of w_p sigma(v_p) v_p (neutrals at rest).
double neutral_ionization_rate(ColliderWithNeutrals& collider, const ParticlesArray& beam) {
    double rate = 0;
    for (int cell = 0; cell < beam.size(); ++cell)
        for (const Particle& p : beam(cell))
            rate += collider.total_collision_frequency(p.velocity, Vector3R::Zero(), beam.mass(), beam.mpw());
    return rate;
}

void run_case(const IonizationCase& c, int steps, double rate_dt, int diag_interval) {
    nlohmann::json domain_config = {{"Dx", 1.0},      {"Dy", 1.0},      {"Dz", 1.0},
                                    {"NumCellsX", 1}, {"NumCellsY", 1}, {"NumCellsZ", 1}};
    const Domain domain = domain_from_config(domain_config);
    const double beam_mass = c.electron_beam ? 1.0 : ion_mass;
    const double beam_energy = (c.electron_beam ? 1.0 : 50.0) / 511.0;   // keV in m_e c^2
    const double speed = std::sqrt(2 * beam_energy / beam_mass);

    Species species;
    // Weight = Density / NumPartPerCell; `count` macroparticles carry `density`.
    species.emplace("Beam", std::make_unique<ParticlesArray>(species_config("Beam", c.electron_beam ? -1.0 : 1.0,
                                                                            beam_mass, c.beam.density, c.beam.count),
                                                             domain));
    species.emplace("Neutrals",
                    std::make_unique<ParticlesArray>(
                        species_config("Neutrals", 0.0, ion_mass, c.neutrals.density, c.neutrals.count), domain));
    species.emplace("Electrons",
                    std::make_unique<ParticlesArray>(
                        species_config("Electrons", -1.0, 1.0, c.electrons.density, c.electrons.count), domain));
    species.emplace("Ions", std::make_unique<ParticlesArray>(
                                species_config("Ions", 1.0, ion_mass, c.ions.density, c.ions.count), domain));
    std::mt19937 engine(2026);
    std::uniform_real_distribution<double> position(0.05, 0.95);
    std::normal_distribution<double> direction(0, 1);
    for (int i = 0; i < c.beam.count; ++i) {
        const Vector3R dir = Vector3R(direction(engine), direction(engine), direction(engine)).normalized();
        species.at("Beam")->add_particle(
            Particle(Vector3R(position(engine), position(engine), position(engine)), speed * dir));
    }
    for (int i = 0; i < c.neutrals.count; ++i)
        species.at("Neutrals")
            ->add_particle(Particle(Vector3R(position(engine), position(engine), position(engine)), Vector3R::Zero()));

    NeutralCollisionConfig config;
    config.charged_species = "Beam";
    config.neutral_species = "Neutrals";
    config.electron_product = "Electrons";
    config.ion_product = "Ions";
    config.electron_ionization = c.electron_beam;
    config.proton_ionization = !c.electron_beam;
    config.proton_charge_exchange = false;
    NeutralCollisionOperator collisions(config, n0);

    CollisionProcessOptions options;
    options.electron_ionization = config.electron_ionization;
    options.proton_ionization = config.proton_ionization;
    options.proton_charge_exchange = false;
    ColliderWithNeutrals rate_model(n0, CollisionScheme::PHYSICAL_ONLY, options);
    const double initial_rate = neutral_ionization_rate(rate_model, *species.at("Beam"));
    const double dt = rate_dt / initial_rate;

    std::ofstream out(c.name + ".csv");
    out << std::setprecision(17) << "step,time,rate,neutrals,electrons,ions,electron_macros,ion_macros\n";
    for (int step = 0; step <= steps; ++step) {
        if (step % diag_interval == 0)
            out << step << "," << step * dt << "," << neutral_ionization_rate(rate_model, *species.at("Beam")) << ","
                << physical(*species.at("Neutrals")) << "," << physical(*species.at("Electrons")) << ","
                << physical(*species.at("Ions")) << "," << macroparticles(*species.at("Electrons")) << ","
                << macroparticles(*species.at("Ions")) << "\n";
        if (step < steps)
            collisions.apply(species, dt);
    }
    nlohmann::json params = {{"case", c.name},
                             {"electron_beam", c.electron_beam},
                             {"dt", dt},
                             {"steps", steps},
                             {"diag_interval", diag_interval},
                             {"weights",
                              {{"beam", species.at("Beam")->mpw()},
                               {"neutrals", species.at("Neutrals")->mpw()},
                               {"electrons", species.at("Electrons")->mpw()},
                               {"ions", species.at("Ions")->mpw()}}},
                             {"neutral_density", c.neutrals.density},
                             {"beam_density", c.beam.density}};
    std::ofstream(c.name + "_params.json") << params.dump(2);
    std::cout << c.name << ": weights beam/neutral " << species.at("Beam")->mpw() / species.at("Neutrals")->mpw()
              << ", electron/neutral " << species.at("Electrons")->mpw() / species.at("Neutrals")->mpw()
              << ", ion/neutral " << species.at("Ions")->mpw() / species.at("Neutrals")->mpw() << ", "
              << macroparticles(*species.at("Neutrals")) << " neutral macroparticles left\n";
}

// Beam of density 1 whose weight is `ratio` times the neutral weight w.
SpeciesWeight beam_weight(double ratio, double w) {
    return {1.0, static_cast<int>(std::lround(1.0 / (ratio * w)))};
}

std::vector<IonizationCase> cases(int n) {
    // Neutral density 0.1, beam density 1 (ten times denser). Product weights are set through Density/ppc.
    const SpeciesWeight neutrals{0.1, n};
    const double w = 0.1 / n;            // neutral weight
    const SpeciesWeight product{w, 1};   // weight w
    const SpeciesWeight double_product{2 * w, 1};
    return {
        {"proton_equal", false, beam_weight(1, w), neutrals, product, product},
        {"proton_beam_heavier_x4", false, beam_weight(4, w), neutrals, product, product},
        {"proton_beam_lighter_x4", false, beam_weight(0.25, w), neutrals, product, product},
        {"proton_dilute_neutrals", false, beam_weight(1, w), {0.05, n}, {0.05 / n, 1}, {0.05 / n, 1}},
        {"electron_equal", true, beam_weight(1, w), neutrals, product, product},
        {"electron_heavy_products", true, beam_weight(1, w), neutrals, double_product, double_product},
        {"electron_unpaired_products", true, beam_weight(1, w), neutrals, product, double_product},
    };
}
}   // namespace

int main(int argc, char** argv) {
    const bool quick = argc == 2 && std::string(argv[1]) == "--quick";
    if (argc != 1 && !quick) {
        std::cerr << "Usage: neutral_ionization [--quick]\n";
        return 1;
    }
    // nu dt = 0.01 at the start; 200 steps decay the neutrals by about e^-2.
    for (const IonizationCase& c : cases(quick ? 4000 : 16000)) run_case(c, 200, 0.01, 5);
    return 0;
}
