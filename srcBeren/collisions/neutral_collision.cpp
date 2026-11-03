#include "collisions/neutral_collision.h"

#include <omp.h>

#include <cstdint>
#include <iostream>

#include "collisions/neutral/collisions_with_neutrals.h"
#include "timer.h"

// Same jump-ahead stride as the legacy collider; each cell borrows a separate slice.
static constexpr int64_t collisionStreamStride = 1 << 16;

NeutralCollisionOperator::NeutralCollisionOperator(const NeutralCollisionConfig& config, double n0)
    : config_(config), n0_(n0), baseRndEng_(77) {
}

std::string NeutralCollisionOperator::info() const {
    return "Neutral(" + config_.charged_species + "-" + config_.neutral_species + ")";
}

void NeutralCollisionOperator::apply(Species& species, const Domain& domain, double dt) {
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

#pragma omp parallel
    {
        timer::flatTimer ompTimer("OMP section");
        ColliderWithNeutrals collider(n0_, scheme, opts);

#pragma omp for schedule(static)
        for (int pk = 0; pk < charged->size(); pk++) {
            int pInCell = (*charged)(pk).size();
            int nInCell = (*neutrals)(pk).size();
            if (pInCell == 0 || nInCell == 0)
                continue;

            LehmerEngine rndEng = baseRndEng_;
            rndEng.discard(static_cast<int64_t>(pk) * collisionStreamStride);

            double n_neutrals = nInCell / (double) neutrals->particles_per_cell();
            std::vector<Particle>& neutrals_data = (*neutrals)(pk);
            int current_neutral_count = nInCell;

            for (int i = 0; i < pInCell && current_neutral_count > 0; i++) {
                std::uniform_int_distribution<> dis(0, current_neutral_count - 1);
                int randomIndex = dis(rndEng);

                Vector3R v1 = (*charged)(pk)[i].velocity;
                Vector3R v2 = neutrals_data[randomIndex].velocity;

                auto [is_collided, ve, vi] = collider.collision_with_neutral(v1, v2, m1, m2, n_neutrals, dt, 1.0 / dt);

                if (is_collided) {
                    Vector3R coord = neutrals_data[randomIndex].coord;
                    Particle pe(coord, ve);
                    Particle pi(coord, vi);

                    electrons->add_particle(pe);
                    ions->add_particle(pi);

                    (*charged)(pk)[i].velocity = v1;
                    std::swap(neutrals_data[randomIndex], neutrals_data[current_neutral_count - 1]);
                    current_neutral_count--;
                } else {
                    (*charged)(pk)[i].velocity = v1;
                    neutrals_data[randomIndex].velocity = v2;
                }
            }

            if (current_neutral_count < static_cast<int>(neutrals_data.size())) {
                neutrals_data.resize(current_neutral_count);
            }
        }
    }

    baseRndEng_.discard(charged->size() * collisionStreamStride);

    neutrals->move(dt);
    neutrals->update_cells(domain);
}
