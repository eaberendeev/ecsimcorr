#include <omp.h>

#include <iomanip>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "collisions/coulomb.h"
#include "test_support.h"

struct Moments {
    Vector3R momentum{0, 0, 0};
    double energy = 0;
    double momentum_scale = 0;   // sum of |particle momentum|: the round-off scale of `momentum`
};

Moments moments(const Species& species) {
    Moments result;
    for (const auto& entry : species) {
        const auto& sp = *entry.second;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const Particle& p : sp(cell)) {
                result.momentum += sp.mpw() * sp.mass() * p.velocity;
                result.momentum_scale += sp.mpw() * sp.mass() * p.velocity.norm();
                result.energy += 0.5 * sp.mpw() * sp.mass() * p.velocity.dot(p.velocity);
            }
    }
    return result;
}

Species mixture(const Domain& domain, int n1, int n2, double density = 1, int ion_ppc = 8) {
    Species species;
    for (const std::string name : {"e", "i"}) {
        auto cfg = check::species(name, name == "e" ? -1 : 1, name == "e" ? 1 : 7, density, name == "e" ? 8 : ion_ppc);
        auto sp = std::make_unique<ParticlesArray>(cfg, domain);
        const int n = name == "e" ? n1 : n2;
        for (int cell = 0; cell < 3; ++cell)
            for (int i = 0; i < n; ++i) {
                const double t = i + cell * 0.7 + (name == "e" ? 0.1 : 0.6);
                sp->add_particle(Particle(Vector3R(0.5 + cell, 0.5, 0.5),
                                          Vector3R(0.01 * std::sin(t), 0.02 * std::cos(t), 0.01 * (i - 2))));
            }
        species.emplace(name, std::move(sp));
    }
    return species;
}

std::vector<Vector3R> velocities(const Species& species) {
    std::vector<Vector3R> result;
    for (const std::string name : {"e", "i"}) {
        const auto& sp = *species.at(name);
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const Particle& p : sp(cell)) result.push_back(p.velocity);
    }
    return result;
}

void axial_relative_velocity() {
    Domain domain;
    domain.init(Vector3I(2, 2, 2), Vector3R(1, 1, 1));
    for (double relative_z : {-0.02, 0.02, 0.0}) {
        const std::string label = "Coulomb axial relative velocity u_z=" + std::to_string(relative_z);
        check::scenario(label, [&] {
            Species species;
            const Vector3R position(0.5, 0.5, 0.5);
            // Equal macro weights and unequal masses test physical momentum as well as the relative-speed rotation.
            for (const std::string name : {"e", "i"}) {
                auto sp = std::make_unique<ParticlesArray>(
                    check::species(name, name == "e" ? -1 : 1, name == "e" ? 1 : 7, 1, 1), domain);
                sp->add_particle(Particle(position, Vector3R(0.012, -0.004, name == "e" ? relative_z : 0.0)));
                species.emplace(name, std::move(sp));
            }
            const auto initial = moments(species);
            const auto before = velocities(species);
            CoulombCollisionOperator collisions("e", "i", 1e10, 15.0);
            collisions.apply(species, 0.03);
            const auto after = velocities(species);
            const auto final = moments(species);

            check::near((after[0] - after[1]).norm(), std::abs(relative_z), 1e-15,
                        "elastic scattering changed relative speed");
            check::near(final.energy, initial.energy, 1e-17, "elastic scattering changed total kinetic energy");
            check::vector(final.momentum, initial.momentum, 1e-15, "elastic scattering changed total momentum");
            for (const auto& entry : species) {
                const auto& sp = *entry.second;
                for (int cell = 0; cell < sp.size(); ++cell)
                    for (const Particle& p : sp(cell)) {
                        check::vector(p.coord, position, 0, "collision moved particle position");
                        for (int d = 0; d < 3; ++d)
                            check::require(std::isfinite(p.velocity[d]), "collision produced nonfinite velocity");
                    }
            }
            if (relative_z == 0) {
                for (size_t i = 0; i < after.size(); ++i)
                    check::vector(after[i], before[i], 0, "zero relative speed changed velocities");
            } else {
                check::require((after[0] - before[0]).norm() > 1e-14,
                               "positive collision rate did not scatter the axial pair");
            }
        });
    }
}

void conservation_and_random_streams() {
    Domain domain;
    domain.init(Vector3I(4, 2, 2), Vector3R(1, 1, 1));
    // Ion weight 1/8, then 1/2 and 1/32 of the electron weight: the moment correction must keep both exact.
    for (int ion_ppc : {8, 2, 32})
        for (int n : {0, 1, 2, 3, 5, 8}) {
            auto species = mixture(domain, n, 7, 1, ion_ppc);
            const auto initial = moments(species);
            CoulombCollisionOperator same("e", "e", 1e10), different("e", "i", 1e10);
            for (int step = 0; step < 4; ++step) {
                same.apply(species, 0.03);
                different.apply(species, 0.03);
            }
            const auto final = moments(species);
            std::cout << "Coulomb conservation: ion ppc=" << ion_ppc << " particles per e-cell=" << n
                      << " |delta E|=" << std::abs(final.energy - initial.energy) << " (limit 1e-14)"
                      << " |delta P|=" << (final.momentum - initial.momentum).norm() << " (component limit 1e-13)\n";
            check::near(final.energy, initial.energy, 1e-14, "Coulomb energy conservation: odd/even/unequal counts");
            check::vector(final.momentum, initial.momentum, 1e-13, "Coulomb momentum conservation");
        }
    auto one = mixture(domain, 5, 3), four = mixture(domain, 5, 3);
    const auto before = velocities(one);
    CoulombCollisionOperator c1("e", "i", 1e10), c4("e", "i", 1e10);
    omp_set_num_threads(1);
    c1.apply(one, 0.03);
    omp_set_num_threads(4);
    c4.apply(four, 0.03);
    auto a = velocities(one), b = velocities(four);
    bool changed = false;
    for (size_t i = 0; i < a.size(); ++i) changed |= (a[i] - before[i]).norm() > 1e-14;
    check::require(changed, "positive collision rate must actually scatter particles");
    for (size_t i = 0; i < a.size(); ++i) check::vector(a[i], b[i], 0, "Coulomb RNG depends on thread count");
    // Density and dt enter only through n*dt: same stream must give the same collision.
    auto base = mixture(domain, 5, 3, 1), dense = mixture(domain, 5, 3, 2);
    CoulombCollisionOperator cb("e", "i", 1e10), cd("e", "i", 1e10);
    cb.apply(base, 0.03);
    cd.apply(dense, 0.015);
    a = velocities(base);
    b = velocities(dense);
    for (size_t i = 0; i < a.size(); ++i) check::vector(a[i], b[i], 1e-15, "Coulomb rate ignores species density");
    cb.apply(base, 0);
    const auto after_zero = velocities(base);
    for (size_t i = 0; i < a.size(); ++i) check::vector(a[i], after_zero[i], 0, "zero timestep changed velocities");
}

// Cold electrons, hot electrons at weights 1/10 and 1/100 of them, and ions; every intra- and interspecies
// channel, including electron populations of equal mass and unequal weight.
void electron_populations_of_unequal_weight_conserve_moments() {
    Domain domain;
    domain.init(Vector3I(3, 1, 1), Vector3R(1, 1, 1));
    const std::vector<std::string> names{"cold", "hot10", "hot100", "ions"};
    const int ppc[] = {4, 40, 400, 4};
    const double mass[] = {1, 1, 1, 7};
    const double thermal[] = {0.002, 0.02, 0.03, 0.001};
    Species species;
    std::mt19937 engine(3);
    std::normal_distribution<double> normal(0, 1);
    for (size_t k = 0; k < names.size(); ++k) {
        auto sp =
            std::make_unique<ParticlesArray>(check::species(names[k], k == 3 ? 1 : -1, mass[k], 1, ppc[k]), domain);
        for (int cell = 0; cell < 3; ++cell)
            for (int i = 0; i < 2 + 3 * cell + static_cast<int>(k); ++i)
                sp->add_particle(Particle(Vector3R(0.5 + cell, 0.5, 0.5),
                                          thermal[k] * Vector3R(normal(engine), normal(engine), normal(engine))));
        species.emplace(names[k], std::move(sp));
    }
    std::vector<std::unique_ptr<CoulombCollisionOperator>> channels;
    for (size_t a = 0; a < names.size(); ++a)
        for (size_t b = a; b < names.size(); ++b)
            channels.push_back(std::make_unique<CoulombCollisionOperator>(names[a], names[b], 1e13));
    const Moments initial = moments(species);
    for (int step = 0; step < 20; ++step)
        for (const std::unique_ptr<CoulombCollisionOperator>& channel : channels) channel->apply(species, 0.5);
    const Moments final = moments(species);
    std::cout << "Coulomb electron populations 1 : 1/10 : 1/100: |delta E|/E="
              << std::abs(final.energy / initial.energy - 1)
              << " |delta P|/sum|p|=" << (final.momentum - initial.momentum).norm() / initial.momentum_scale << "\n";
    check::near(final.energy, initial.energy, 1e-13 * initial.energy, "electron populations energy");
    check::vector(final.momentum, initial.momentum, 1e-14 * initial.momentum_scale, "electron populations momentum");
}

// Temperature of species `name` (mass times mean squared velocity per degree of freedom).
double temperature(const Species& species, const std::string& name) {
    const ParticlesArray& sp = *species.at(name);
    double sum = 0;
    int count = 0;
    for (int cell = 0; cell < sp.size(); ++cell)
        for (const Particle& p : sp(cell)) {
            sum += p.velocity.dot(p.velocity);
            ++count;
        }
    return sp.mass() * sum / (3.0 * count);
}

// Hot electrons and cold ions in one cell: the ion heating after a few steps must not depend on the
// ion weight beyond sampling noise (rejection scattering is correct on average).
void unequal_weights_relax_like_equal_weights() {
    Domain domain;
    domain.init(Vector3I(1, 1, 1), Vector3R(1, 1, 1));
    std::mt19937 engine(7);
    std::normal_distribution<double> normal(0, 1);
    double heating[2] = {0, 0};
    double electron_temperature = 0, ion_temperature = 0;
    const int ion_ppc[2] = {8000, 2000};   // equal weights, then ion weight 4x the electron weight
    for (int run = 0; run < 2; ++run) {
        Species species;
        for (const std::string name : {"e", "i"}) {
            const bool electron = name == "e";
            const int ppc = electron ? 8000 : ion_ppc[run];
            auto sp = std::make_unique<ParticlesArray>(
                check::species(name, electron ? -1 : 1, electron ? 1 : 4, 1000, ppc), domain);
            const double thermal = electron ? 0.02 : 0.002;
            for (int i = 0; i < ppc; ++i)
                sp->add_particle(Particle(Vector3R(0.5, 0.5, 0.5),
                                          thermal * Vector3R(normal(engine), normal(engine), normal(engine))));
            species.emplace(name, std::move(sp));
        }
        const double t0 = temperature(species, "i");
        electron_temperature = temperature(species, "e");
        ion_temperature = t0;
        const auto initial = moments(species);
        CoulombCollisionOperator collisions("e", "i", 1e13);
        for (int step = 0; step < 400; ++step) collisions.apply(species, 0.5);
        const auto final = moments(species);
        check::near(final.energy, initial.energy, 1e-13 * initial.energy, "weighted relaxation energy");
        std::cout << "weighted relaxation: run " << run
                  << " |delta P| / sum |p| = " << (final.momentum - initial.momentum).norm() / initial.momentum_scale
                  << "\n";
        check::vector(final.momentum, initial.momentum, 1e-14 * initial.momentum_scale, "weighted relaxation momentum");
        heating[run] = temperature(species, "i") - t0;
    }
    std::cout << "Coulomb ion heating: equal weights " << heating[0] << ", ion weight x4 " << heating[1]
              << " (initial gap " << electron_temperature - ion_temperature << ")\n";
    check::require(heating[0] > 0, "electron-ion relaxation must heat the cold ions");
    // One run scatters by about 3% (eight seeds: mean -1.6%); 10% is about three standard deviations and
    // still catches a wrong weight factor such as a missing w_min / w_max.
    check::near(heating[1], heating[0], 0.1 * heating[0], "unequal weights change the relaxation rate");
}

int main() {
    return check::run([] {
        std::cout << std::scientific << std::setprecision(8);
        axial_relative_velocity();
        conservation_and_random_streams();
        unequal_weights_relax_like_equal_weights();
        electron_populations_of_unequal_weight_conserve_moments();
        std::cout << "PASS: Coulomb conservation, density scaling and random streams\n";
    });
}
