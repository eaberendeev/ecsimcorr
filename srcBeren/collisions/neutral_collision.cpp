#include "collisions/neutral_collision.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "collisions/neutral/collisions_with_neutrals.h"
#include "runtime/profiler.h"

// Each cell borrows a separate slice of the process random stream.
static constexpr int64_t collisionStreamStride = 1 << 16;

namespace {
// Macroparticles of weight `weight` carrying `events` physical events: the integer part, plus one more with the
// probability of the remainder, as for secondary emission. An integer ratio draws no random number.
int macro_count(double events, double weight, std::uniform_real_distribution<double>& uniform, LehmerEngine& engine) {
    const double expected = events / weight;
    const double whole = std::floor(expected);
    if (expected == whole)
        return static_cast<int>(whole);
    return static_cast<int>(whole) + (uniform(engine) < expected - whole ? 1 : 0);
}

// Swap the neutral to the end of the live range and shorten it.
void remove_neutral(std::vector<Particle>& neutrals, int index, int& live) {
    std::swap(neutrals[index], neutrals[live - 1]);
    --live;
}
}   // namespace

NeutralCollisionOperator::NeutralCollisionOperator(const NeutralCollisionConfig& config, double n0)
    : config_(config), n0_(n0), baseRndEng_(77) {
}

std::string NeutralCollisionOperator::info() const {
    return "Neutral(" + config_.charged_species + "-" + config_.neutral_species + ")";
}

void NeutralCollisionOperator::apply(Species& species, double dt) {
    RECORD_TIMER;

    ParticlesArray* charged = find_species(species, config_.charged_species);
    ParticlesArray* neutrals = find_species(species, config_.neutral_species);
    ParticlesArray* electrons = find_species(species, config_.electron_product);
    ParticlesArray* ions = find_species(species, config_.ion_product);

    if (!charged || !neutrals || !electrons || !ions) {
        std::cerr << "NeutralCollision: missing species\n";
        return;
    }

    CollisionScheme scheme = CollisionScheme::PHYSICAL_ONLY;
    if (config_.scheme == "null_collision") {
        scheme = CollisionScheme::NULL_COLLISION;
    }

    CollisionProcessOptions opts;
    opts.electron_ionization = config_.electron_ionization;
    opts.proton_ionization = config_.proton_ionization;
    opts.proton_charge_exchange = config_.proton_charge_exchange;

    const double m1 = charged->mass();
    const double m2 = neutrals->mass();
    // A reaction of one projectile macroparticle stands for w_charged physical reactions; neutrals, electrons and
    // ions of other weights take part with w_charged / w_species macroparticles. Equal product weights create
    // electrons and ions in pairs, which keeps the charge exact.
    const double w_charged = charged->mpw();
    const double w_neutral = neutrals->mpw();
    const double w_electron = electrons->mpw();
    const double w_ion = ions->mpw();
    const bool paired_products = w_electron == w_ion;

#pragma omp parallel
    {
        timer::FlatTimer ompTimer("OMP section");
        ColliderWithNeutrals collider(n0_, scheme, opts);

#pragma omp for schedule(static)
        for (int pk = 0; pk < charged->size(); pk++) {
            int pInCell = (*charged)(pk).size();
            int nInCell = (*neutrals)(pk).size();
            if (pInCell == 0 || nInCell == 0)
                continue;

            LehmerEngine rndEng = baseRndEng_;
            rndEng.discard(static_cast<int64_t>(pk) * collisionStreamStride);

            const double n_neutrals = nInCell * w_neutral;   // physical neutral density in units of n0
            std::vector<Particle>& neutrals_data = (*neutrals)(pk);
            int current_neutral_count = nInCell;
            std::uniform_real_distribution<double> uniform(0.0, 1.0);

            for (int i = 0; i < pInCell && current_neutral_count > 0; i++) {
                std::uniform_int_distribution<> dis(0, current_neutral_count - 1);
                int randomIndex = dis(rndEng);

                Vector3R v1 = (*charged)(pk)[i].velocity;
                Vector3R v2 = neutrals_data[randomIndex].velocity;
                const Vector3R neutral_velocity = v2;

                const auto [is_collided, ve, vi] =
                    collider.collision_with_neutral(v1, v2, m1, m2, n_neutrals, dt, 1.0 / dt);
                (*charged)(pk)[i].velocity = v1;
                const bool exchanged = !is_collided && v2 != neutral_velocity;
                if (!is_collided && !exchanged)
                    continue;

                const int neutrals_hit =
                    std::min(current_neutral_count, macro_count(w_charged, w_neutral, uniform, rndEng));
                if (exchanged) {
                    // Charge exchange: distinct hit neutrals take the former projectile velocity. Extra hits are
                    // drawn without repetition by moving the chosen neutrals behind the drawing range.
                    if (neutrals_hit > 0)
                        neutrals_data[randomIndex].velocity = v2;
                    int range = current_neutral_count;
                    if (neutrals_hit > 1)
                        std::swap(neutrals_data[randomIndex], neutrals_data[--range]);
                    for (int extra = 1; extra < neutrals_hit; ++extra) {
                        std::swap(neutrals_data[std::uniform_int_distribution<>(0, range - 1)(rndEng)],
                                  neutrals_data[range - 1]);
                        neutrals_data[--range].velocity = v2;
                    }
                    continue;
                }
                const Vector3R coord = neutrals_data[randomIndex].coord;
                const int new_electrons = macro_count(w_charged, w_electron, uniform, rndEng);
                const int new_ions = paired_products ? new_electrons : macro_count(w_charged, w_ion, uniform, rndEng);
                for (int k = 0; k < new_electrons; ++k) electrons->add_particle(Particle(coord, ve));
                for (int k = 0; k < new_ions; ++k) ions->add_particle(Particle(coord, vi));
                if (neutrals_hit > 0)
                    remove_neutral(neutrals_data, randomIndex, current_neutral_count);
                for (int extra = 1; extra < neutrals_hit; ++extra)
                    remove_neutral(neutrals_data, std::uniform_int_distribution<>(0, current_neutral_count - 1)(rndEng),
                                   current_neutral_count);
            }

            if (current_neutral_count < static_cast<int>(neutrals_data.size())) {
                neutrals_data.resize(current_neutral_count);
            }
        }
    }

    baseRndEng_.discard(charged->size() * collisionStreamStride);
}
