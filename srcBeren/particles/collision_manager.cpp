#include "collision_manager.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>

#include "collisions_with_neutrals.h"
#include "vector3.h"

static inline double get_reduced_mass(double m1, double m2) {
    return m1 * m2 / (m1 + m2);
}

// ============================================================================
// CoulombCollisionOperator
// ============================================================================

CoulombCollisionOperator::CoulombCollisionOperator(const std::string &species1, const std::string &species2, double n0,
                                                   double coulomb_log)
    : species1_name_(species1),
      species2_name_(species2),
      is_same_type_(species1 == species2),
      n0_(n0),
      coulomb_log_(coulomb_log) {
    gen_.SetRandSeed(42);
}

std::string CoulombCollisionOperator::info() const {
    if (is_same_type_) {
        return "Coulomb(" + species1_name_ + "-" + species2_name_ + ", lnL=" + std::to_string(coulomb_log_) + ")";
    }
    return "Coulomb(" + species1_name_ + "-" + species2_name_ + ", lnL=" + std::to_string(coulomb_log_) + ")";
}

double CoulombCollisionOperator::get_variance(double u, double q1, double q2, double n, double m, double dt) {
    assert(u > 0);
    const double wp3 = pow(SGS::get_plasma_freq(n0_), 3);
    constexpr double c3 = SGS::c * SGS::c * SGS::c;
    return (wp3 / (c3 * n0_)) * (coulomb_log_ * q1 * q1 * q2 * q2 * n * dt) / (8.0 * M_PI * m * m * u * u * u);
}

void CoulombCollisionOperator::bin_collide(Vector3R &v1, Vector3R &v2, double q1, double q2, double n1, double n2,
                                           double m1, double m2, double dt, double variance_factor,
                                           ThreadRandomGenerator &rng) {
    const double n = std::min(n1, n2);
    const double m = get_reduced_mass(m1, m2);
    const Vector3R u = v1 - v2;
    const double modu = u.norm();
    if (modu < 1.e-30)
        return;

    const double variance = variance_factor * get_variance(modu, q1, q2, n, m, dt);
    const double sigma = (variance < 1.0) ? rng.Gauss(sqrt(variance)) : M_PI * rng.Uniform01();
    const double phi = 2.0 * M_PI * rng.Uniform01();
    const double cosp = cos(phi);
    const double sinp = sin(phi);
    const double sint = 2.0 * sigma / (1.0 + sigma * sigma);
    const double cost = 1.0 - 2.0 * sigma * sigma / (1.0 + sigma * sigma);
    const double up = sqrt(u.x() * u.x() + u.y() * u.y());

    Vector3R du;
    // Two branches for the Takizuka-Abe rotation:
    //
    // Branch 1 (up ≈ 0): u is nearly parallel to z — rotate in XY plane.
    //   du_xy = modu * sin(θ) * (cos φ, sin φ),  du_z = -sign(u_z) * modu * (1 - cos θ)
    //   Result: |u'| = |u|, u' is rotated by θ isotropically around z.
    //
    // Branch 2 (general 3D): standard rotation of u by θ around an axis
    //   perpendicular to u at azimuth φ.  When up → 0, branch 2 approaches:
    //     dux ≈ modu * sin θ * cos(α + φ),  duy ≈ modu * sin θ * sin(α + φ)
    //   where α is the azimuth of the vanishing perpendicular component.
    //   Since φ ~ Uniform[0, 2π), the phase shift α is irrelevant —
    //   both branches give the same distribution.
    if (up < 1.e-16 * modu) {
        double sign_z = (u.z() >= 0.0) ? 1.0 : -1.0;
        du = Vector3R(modu * sint * cosp, modu * sint * sinp, -sign_z * modu * (1.0 - cost));
    } else {
        du = Vector3R((u.x() / up) * u.z() * sint * cosp - (u.y() / up) * modu * sint * sinp - u.x() * (1.0 - cost),
                      (u.y() / up) * u.z() * sint * cosp + (u.x() / up) * modu * sint * sinp - u.y() * (1.0 - cost),
                      -up * sint * cosp - u.z() * (1.0 - cost));
    }

    const double mu_over_m1 = m / m1;
    const double mu_over_m2 = m / m2;
    v1 += mu_over_m1 * du;
    v2 -= mu_over_m2 * du;
}

void CoulombCollisionOperator::collide_same_type(ParticlesArray &sp, double dt) {
    const double q = sp.charge;
    const double m1 = sp.mass();
    const int num_cells = sp.size();

#pragma omp parallel
    {
        std::vector<int> indices;

#pragma omp for schedule(static)
        for (int pk = 0; pk < num_cells; pk++) {
            const int N = sp.particlesData(pk).size();
            if (N < 2)
                continue;

            indices.resize(N);
            std::iota(indices.begin(), indices.end(), 0);

            auto &rng_engine = gen_.gen();
            const bool is_odd = (N % 2 == 1);
            const double n_density = N / (double) sp.NumPartPerCell;

            if (is_odd) {
                // Case 1b: first 3 particles form 3 pairs with variance_factor=0.5
                // Lazy Fisher-Yates for first 3 elements
                for (int k = 0; k < 3 && k < N; k++) {
                    std::uniform_int_distribution<int> dist(k, N - 1);
                    std::swap(indices[k], indices[dist(rng_engine)]);
                }

                auto do_collide = [&](int i1, int i2, double vf) {
                    bin_collide(sp.particlesData(pk)[i1].velocity, sp.particlesData(pk)[i2].velocity, q, q, n_density,
                                n_density, m1, m1, dt, vf, gen_);
                };

                do_collide(indices[0], indices[1], 0.5);
                do_collide(indices[1], indices[2], 0.5);
                do_collide(indices[0], indices[2], 0.5);

                // Remaining pairs: lazy Fisher-Yates from index 3 onward
                for (int k = 3; k + 1 < N; k += 2) {
                    std::uniform_int_distribution<int> dist1(k, N - 1);
                    std::swap(indices[k], indices[dist1(rng_engine)]);
                    std::uniform_int_distribution<int> dist2(k + 1, N - 1);
                    std::swap(indices[k + 1], indices[dist2(rng_engine)]);
                    do_collide(indices[k], indices[k + 1], 1.0);
                }
            } else {
                // Case 1a: all pairs with variance_factor=1.0
                // Lazy Fisher-Yates: shuffle only as needed
                for (int k = 0; k + 1 < N; k += 2) {
                    std::uniform_int_distribution<int> dist1(k, N - 1);
                    std::swap(indices[k], indices[dist1(rng_engine)]);
                    std::uniform_int_distribution<int> dist2(k + 1, N - 1);
                    std::swap(indices[k + 1], indices[dist2(rng_engine)]);

                    bin_collide(sp.particlesData(pk)[indices[k]].velocity,
                                sp.particlesData(pk)[indices[k + 1]].velocity, q, q, n_density, n_density, m1, m1, dt,
                                1.0, gen_);
                }
            }
        }
    }
}

void CoulombCollisionOperator::collide_diff_type(ParticlesArray &sp1, ParticlesArray &sp2, double dt) {
    const double q1 = sp1.charge;
    const double m1 = sp1.mass();
    const double q2 = sp2.charge;
    const double m2 = sp2.mass();
    const int num_cells = sp1.size();

#pragma omp parallel
    {
        std::vector<int> indices_large;
        std::vector<int> indices_small;

#pragma omp for schedule(static)
        for (int pk = 0; pk < num_cells; pk++) {
            const int N1 = sp1.particlesData(pk).size();
            const int N2 = sp2.particlesData(pk).size();
            if (N1 == 0 || N2 == 0)
                continue;

            const bool sp1_is_larger = (N1 >= N2);
            const int N_large = sp1_is_larger ? N1 : N2;
            const int N_small = sp1_is_larger ? N2 : N1;

            indices_large.resize(N_large);
            std::iota(indices_large.begin(), indices_large.end(), 0);

            indices_small.resize(N_small);
            std::iota(indices_small.begin(), indices_small.end(), 0);

            auto &rng_engine = gen_.gen();

            // Lazy Fisher-Yates for the large array
            // For the small array, we cycle through deterministically (Takizuka Case 2)
            const int r = N_large % N_small;
            const int quot = N_large / N_small;
            const int first_group_size = r * (quot + 1);

            // Shuffle the large array lazily as we iterate
            // Shuffle the small array once (full, since all get reused)
            for (int k = 0; k < N_small; k++) {
                std::uniform_int_distribution<int> dist(k, N_small - 1);
                std::swap(indices_small[k], indices_small[dist(rng_engine)]);
            }

            const double n1_density = N1 / (double) sp1.NumPartPerCell;
            const double n2_density = N2 / (double) sp2.NumPartPerCell;

            for (int first_ind = 0; first_ind < N_large; first_ind++) {
                // Lazy Fisher-Yates for large array
                std::uniform_int_distribution<int> dist(first_ind, N_large - 1);
                std::swap(indices_large[first_ind], indices_large[dist(rng_engine)]);

                // Determine partner from small array (Takizuka Case 2 mapping)
                int second_ind;
                if (first_ind < first_group_size) {
                    second_ind = first_ind / (quot + 1);
                } else {
                    second_ind = r + (first_ind - first_group_size) / quot;
                }

                int idx1, idx2;
                if (sp1_is_larger) {
                    idx1 = indices_large[first_ind];
                    idx2 = indices_small[second_ind];
                } else {
                    idx1 = indices_small[second_ind];
                    idx2 = indices_large[first_ind];
                }

                bin_collide(sp1.particlesData(pk)[idx1].velocity, sp2.particlesData(pk)[idx2].velocity, q1, q2,
                            n1_density, n2_density, m1, m2, dt, 1.0, gen_);
            }
        }
    }
}

void CoulombCollisionOperator::apply(Species &species, [[maybe_unused]] const Domain &domain, double dt) {
    if (is_same_type_) {
        ParticlesArray *sp = find_species(species, species1_name_);
        if (!sp) {
            std::cerr << "CoulombCollision: species '" << species1_name_ << "' not found\n";
            return;
        }
        collide_same_type(*sp, dt);
    } else {
        ParticlesArray *sp1 = find_species(species, species1_name_);
        ParticlesArray *sp2 = find_species(species, species2_name_);
        if (!sp1 || !sp2) {
            std::cerr << "CoulombCollision: species '" << species1_name_ << "' or '" << species2_name_
                      << "' not found\n";
            return;
        }
        collide_diff_type(*sp1, *sp2, dt);
    }
}

// ============================================================================
// NeutralCollisionOperator
// ============================================================================

NeutralCollisionOperator::NeutralCollisionOperator(const NeutralCollisionConfig &config, double n0)
    : config_(config), n0_(n0) {
    gen_.SetRandSeed(77);
}

std::string NeutralCollisionOperator::info() const {
    return "Neutral(" + config_.charged_species + "-" + config_.neutral_species + ")";
}

void NeutralCollisionOperator::apply(Species &species, const Domain &domain, double dt) {
    ParticlesArray *charged = find_species(species, config_.charged_species);
    ParticlesArray *neutrals = find_species(species, config_.neutral_species);
    ParticlesArray *electrons = find_species(species, config_.electron_product);
    ParticlesArray *ions = find_species(species, config_.ion_product);

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
        ColliderWithNeutrals collider(n0_, scheme, opts);

#pragma omp for schedule(static)
        for (int pk = 0; pk < charged->size(); pk++) {
            int pInCell = charged->particlesData(pk).size();
            int nInCell = neutrals->particlesData(pk).size();
            if (pInCell == 0 || nInCell == 0)
                continue;

            double n1 = pInCell / (double) charged->NumPartPerCell;
            auto &neutrals_data = neutrals->particlesData(pk);
            int current_neutral_count = nInCell;

            for (int i = 0; i < pInCell && current_neutral_count > 0; i++) {
                std::uniform_int_distribution<> dis(0, current_neutral_count - 1);
                int randomIndex = dis(gen_.gen());

                Particle &charged_particle = charged->particlesData(pk)[i];
                Particle &neutral_particle = neutrals_data[randomIndex];

                Vector3R v1 = charged_particle.velocity;
                Vector3R v2 = neutral_particle.velocity;

                auto [is_collided, ve, vi] = collider.collision_with_neutral(v1, v2, m1, m2, n1, dt, 1.0 / dt);

                if (is_collided) {
                    Vector3R coord = neutral_particle.coord;
                    Particle pe(coord, ve);
                    Particle pi(coord, vi);

                    electrons->add_particle(pe);
                    ions->add_particle(pi);

                    charged_particle.velocity = v1;
                    std::swap(neutral_particle, neutrals_data[current_neutral_count - 1]);
                    current_neutral_count--;
                } else {
                    charged_particle.velocity = v1;
                    neutral_particle.velocity = v2;
                }
            }

            if (current_neutral_count < static_cast<int>(neutrals_data.size())) {
                neutrals_data.resize(current_neutral_count);
            }
        }
    }

    neutrals->move(dt);
    neutrals->update_cells(domain);
}

// ============================================================================
// CollisionManager
// ============================================================================

void CollisionManager::init_from_json(const nlohmann::json &config, double n0) {
    operators_.clear();

    if (!config.contains("Collisions")) {
        return;
    }

    const auto &collisions = config["Collisions"];
    if (!collisions.is_array()) {
        std::cerr << "CollisionManager: 'Collisions' must be an array\n";
        return;
    }

    for (const auto &entry : collisions) {
        if (!entry.contains("type")) {
            std::cerr << "CollisionManager: entry missing 'type'\n";
            continue;
        }

        const std::string type = entry["type"].get<std::string>();

        if (type == "coulomb") {
            std::string sp1 = entry.value("species1", "");
            std::string sp2 = entry.value("species2", "");
            double lnL = entry.value("coulomb_log", 15.0);
            if (sp1.empty() || sp2.empty()) {
                std::cerr << "CollisionManager: coulomb entry needs species1, species2\n";
                continue;
            }
            operators_.push_back(std::make_unique<CoulombCollisionOperator>(sp1, sp2, n0, lnL));
        } else if (type == "neutral_ionization") {
            NeutralCollisionConfig nc;
            nc.charged_species = entry.value("charged_species", "");
            nc.neutral_species = entry.value("neutral_species", "Neutrals");
            nc.electron_product = entry.value("electron_product", "Electrons");
            nc.ion_product = entry.value("ion_product", "Ions");
            nc.scheme = entry.value("scheme", "physical_only");
            nc.electron_ionization = entry.value("electron_ionization", true);
            nc.proton_ionization = entry.value("proton_ionization", true);
            nc.proton_charge_exchange = entry.value("proton_charge_exchange", true);
            if (nc.charged_species.empty()) {
                std::cerr << "CollisionManager: neutral_ionization needs charged_species\n";
                continue;
            }
            operators_.push_back(std::make_unique<NeutralCollisionOperator>(nc, n0));
        } else {
            std::cerr << "CollisionManager: unknown collision type '" << type << "'\n";
        }
    }

    if (operators_.empty())
        return;
    std::cout << "Collisions:\n";
    for (const auto &op : operators_) {
        std::cout << "  " << op->info() << "\n";
    }
}

void CollisionManager::apply(Species &species, const Domain &domain, double dt) {
    for (auto &op : operators_) {
        op->apply(species, domain, dt);
    }
}
