#include <omp.h>

#include "collision_manager.h"
#include "test_support.h"

struct Moments {
    Vector3R momentum{0, 0, 0};
    double energy = 0;
};
Moments moments(const Species& species) {
    Moments result;
    for (const auto& entry : species) {
        const auto& sp = *entry.second;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const auto& p : sp.particlesData(cell)) {
                result.momentum += sp.mpw() * sp.mass() * p.velocity;
                result.energy += 0.5 * sp.mpw() * sp.mass() * p.velocity.dot(p.velocity);
            }
    }
    return result;
}
Species mixture(const Domain& domain, int n1, int n2, double density = 1) {
    Species species;
    for (const std::string name : {"e", "i"}) {
        auto cfg = check::species(name, name == "e" ? -1 : 1, name == "e" ? 1 : 7, density, 8);
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
            for (const auto& p : sp.particlesData(cell)) result.push_back(p.velocity);
    }
    return result;
}
void test_collisions() {
    Domain domain;
    domain.init(Vector3I(4, 2, 2), Vector3R(1, 1, 1));
    for (int n : {0, 1, 2, 3, 5, 8}) {
        auto species = mixture(domain, n, 7);
        const auto initial = moments(species);
        CoulombCollisionOperator same("e", "e", 1e10), different("e", "i", 1e10);
        for (int step = 0; step < 4; ++step) {
            same.apply(species, domain, 0.03);
            different.apply(species, domain, 0.03);
        }
        const auto final = moments(species);
        check::near(final.energy, initial.energy, 1e-14, "Coulomb energy conservation: odd/even/unequal counts");
        check::vector(final.momentum, initial.momentum, 1e-13, "Coulomb momentum conservation");
    }
    auto one = mixture(domain, 5, 3), four = mixture(domain, 5, 3);
    const auto before = velocities(one);
    CoulombCollisionOperator c1("e", "i", 1e10), c4("e", "i", 1e10);
    omp_set_num_threads(1);
    c1.apply(one, domain, 0.03);
    omp_set_num_threads(4);
    c4.apply(four, domain, 0.03);
    auto a = velocities(one), b = velocities(four);
    bool changed = false;
    for (size_t i = 0; i < a.size(); ++i) changed |= (a[i] - before[i]).norm() > 1e-14;
    check::require(changed, "positive collision rate must actually scatter particles");
    for (size_t i = 0; i < a.size(); ++i) check::vector(a[i], b[i], 0, "Coulomb RNG depends on thread count");
    // Density and dt enter only through n*dt: same stream must give the same collision.
    auto base = mixture(domain, 5, 3, 1), dense = mixture(domain, 5, 3, 2);
    CoulombCollisionOperator cb("e", "i", 1e10), cd("e", "i", 1e10);
    cb.apply(base, domain, 0.03);
    cd.apply(dense, domain, 0.015);
    a = velocities(base);
    b = velocities(dense);
    for (size_t i = 0; i < a.size(); ++i) check::vector(a[i], b[i], 1e-15, "Coulomb rate ignores species density");
    cb.apply(base, domain, 0);
    const auto after_zero = velocities(base);
    for (size_t i = 0; i < a.size(); ++i) check::vector(a[i], after_zero[i], 0, "zero timestep changed velocities");
    auto unequal = mixture(domain, 5, 3);
    unequal.at("i") = std::make_unique<ParticlesArray>(check::species("i", 1, 7, 2, 8), domain);
    bool rejected = false;
    try {
        cb.apply(unequal, domain, 0.03);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "unsupported unequal weights silently violate physical momentum conservation");
}
int main() {
    return check::run(test_collisions);
}
