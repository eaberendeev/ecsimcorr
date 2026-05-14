// #include "Particles.h"

#include "collision.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "collisions_with_neutrals.h"
#include "vector3.h"

double BinaryCollider::get_variance_coll(double u, double q1, double q2, double n, double m, double dt) {
    const double lk = 15;
    return (pow(SGS::get_plasma_freq(n0), 3) / (SGS::c * SGS::c * SGS::c * n0)) * (lk * q1 * q1 * q2 * q2 * n * dt) /
           (8 * M_PI * m * m * u * u * u);
}

std::pair<int, int> BinaryCollisionSameType::get_pair() {
    // Check correct particles size
    if (!canCollide()) {
        std::cout << "Event is not exist\n";
        return std::pair<int, int>(-1, -1);
    }

    event++;
    if (isOddEvents) {
        if (event == 0) {
            return std::pair<int, int>(v[0], v[1]);
        } else if (event == 1) {
            return std::pair<int, int>(v[1], v[2]);
        } else if (event == 2) {
            return std::pair<int, int>(v[0], v[2]);
        }
        return std::pair<int, int>(v[2 * event - 3], v[2 * event - 2]);
    }

    return std::pair<int, int>(v[2 * event], v[2 * event + 1]);
}

std::pair<int, int> BinaryCollisionDiffType::get_pair() {
    if (!canCollide()) {
        std::cout << "Event is not exist\n";
        return std::pair<int, int>(-1, -1);
    }
    event++;

    if (v1.size() == v2.size()) {
        return std::pair<int, int>(v1[event], v2[event]);
    }

    const int r = v1.size() % v2.size();
    const int q = v1.size() / v2.size();
    const int firstGroup = r * (q + 1);   // number of particles in first group in sort1

    int secondInd;

    if (firstInd < firstGroup) {
        secondInd = firstInd / (q + 1);
    } else {
        secondInd = r + (firstInd - firstGroup) / q;
    }
    int ind1 = firstInd;
    int ind2 = secondInd;

    firstInd++;

    // swap sort of particles
    if (leadingSort == 1) {
        return std::pair<int, int>(v2[ind2], v1[ind1]);
    }

    return std::pair<int, int>(v1[ind1], v2[ind2]);
}

void BinaryCollider::bin_collide(Vector3R &v1, Vector3R &v2, double q1, double q2, double n1, double n2, double m1,
                                 double m2, double dt, double variance_factor) {
    const double n = std::min(n1, n2);
    const double m = get_center_mass(m1, m2);
    const Vector3R u = v1 - v2;
    const double modu = u.norm();
    const double variance = variance_factor * get_variance_coll(modu, q1, q2, n, m, dt);

    const double sigma = (variance < 1) ? gen.Gauss(sqrt(variance)) : M_PI * gen.Uniform01();

    const double phi = 2 * M_PI * gen.Uniform01();
    const double cosp = cos(phi);
    const double sinp = sin(phi);
    const double sint = 2 * sigma / (1 + sigma * sigma);
    const double cost = 1 - 2 * sigma * sigma / (1 + sigma * sigma);
    const double up = sqrt(u.x() * u.x() + u.y() * u.y());
    double dux, duy, duz;
    if (up < 1.e-16) {
        dux = modu * sint * cosp;
        duy = modu * sint * sinp;
        duz = -modu * (1 - cost);
    } else {
        dux = (u.x() / up) * u.z() * sint * cosp - (u.y() / up) * modu * sint * sinp - u.x() * (1 - cost);
        duy = (u.y() / up) * u.z() * sint * cosp + (u.x() / up) * modu * sint * sinp - u.y() * (1 - cost);
        duz = -up * sint * cosp - u.z() * (1 - cost);
    }
    const Vector3R du = Vector3R(dux, duy, duz);
    v1 += (m / m1) * du;
    v2 -= (m / m2) * du;
}

void BinaryCollider::collide_same_sort_binary(Species &species, const double dt) {
    for (auto &kv : species) {
        auto &sp = *kv.second;
        const double q = sp.charge;
        const double m1 = sp.mass();
#pragma omp parallel for schedule(dynamic, 32)
        for (auto pk = 0; pk < sp.size(); pk++) {
            BinaryCollisionSameType collider(sp.particlesData(pk).size(), gen.gen());
            while (collider.canCollide()) {
                auto pair = collider.get_pair();
                Vector3R v1 = sp.particlesData(pk)[pair.first].velocity;
                Vector3R v2 = sp.particlesData(pk)[pair.second].velocity;
                double n1 = sp.particlesData(pk).size() / (double) sp.NumPartPerCell;
                double variance_factor = collider.get_variance_factor();
                bin_collide(v1, v2, q, q, n1, n1, m1, m1, dt, variance_factor);
                sp.particlesData(pk)[pair.first].velocity = v1;
                sp.particlesData(pk)[pair.second].velocity = v2;
            }
        }
    }
}

void BinaryCollider::collide_ion_electron_binary(Species &species, const double dt) {
    auto *e = find_species(species, "Electrons");
    auto *i = find_species(species, "Ions");
    if (!e || !i) {
        std::cout << "Error: Invalid type of species for collision\n";
        return;
    }

    const double q1 = e->charge;
    const double m1 = e->mass();
    const double q2 = i->charge;
    const double m2 = i->mass();
#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < e->size(); pk++) {
        BinaryCollisionDiffType collider(e->particlesData(pk).size(), i->particlesData(pk).size(), gen.gen());
        while (collider.canCollide()) {
            auto pair = collider.get_pair();
            Vector3R v1 = e->particlesData(pk)[pair.first].velocity;
            Vector3R v2 = i->particlesData(pk)[pair.second].velocity;
            const double variance_factor = 1.;
            double n1 = e->particlesData(pk).size() / (double) e->NumPartPerCell;
            double n2 = i->particlesData(pk).size() / (double) i->NumPartPerCell;
            bin_collide(v1, v2, q1, q2, n1, n2, m1, m2, dt, variance_factor);
            e->particlesData(pk)[pair.first].velocity = v1;
            i->particlesData(pk)[pair.second].velocity = v2;
        }
    }
}

void BinaryColliderWithNeutrals::collide_with_neutrals_binary(Species &species, const Domain &domain, const double dt) {
    for (auto &kv : species) {
        auto &sp = *kv.second;
        if (sp.is_neutral())
            continue;
        collide_with_neutrals_binary_impl(species, kv.first, dt);
    }

    // Move neutrals
    if (auto *neutrals = find_species(species, "Neutrals")) {
        neutrals->move(dt);
        neutrals->update_cells(domain);
    }
}

void BinaryColliderWithNeutrals::collide_with_neutrals_binary_impl(Species &species, const std::string &pType,
                                                                   const double dt) {
    auto *neutrals = find_species(species, "Neutrals");
    auto *electrons = find_species(species, "Electrons");
    auto *ions = find_species(species, "Ions");
    auto *p = find_species(species, pType);

    if (!neutrals || !electrons || !ions || !p) {
        std::cout << "Error: Invalid type of species for collision\n";
        return;
    }

    const double m1 = p->mass();
    const double m2 = neutrals->mass();

#pragma omp parallel
    {
        ColliderWithNeutrals colliderWithNeutrals(n0, scheme, process_opts);

#pragma omp for schedule(dynamic, 32)
        for (auto pk = 0; pk < p->size(); pk++) {
            int pInCell = p->particlesData(pk).size();
            int nInCell = neutrals->particlesData(pk).size();

            if (pInCell == 0 || nInCell == 0)
                continue;

            double n1 = pInCell / (double) p->NumPartPerCell;
            // double n2 = nInCell / (double)
            // species[neutrals_type]->NumPartPerCell;

            // Работаем напрямую с данными нейтралов
            auto &neutrals_data = neutrals->particlesData(pk);
            int current_neutral_count = nInCell;

            for (int i = 0; i < pInCell && current_neutral_count > 0; i++) {
                std::uniform_int_distribution<> dis(0, current_neutral_count - 1);
                int randomIndex = dis(gen.gen());

                Particle &charged_particle = p->particlesData(pk)[i];
                Particle &neutral_particle = neutrals_data[randomIndex];

                Vector3R v1 = charged_particle.velocity;
                Vector3R v2 = neutral_particle.velocity;

                auto [is_collided, ve, vi] =
                    colliderWithNeutrals.collision_with_neutral(v1, v2, m1, m2, n1, dt, 1. / dt);
                if (is_collided) {
                    Vector3R coord = neutral_particle.coord;
                    Particle pe(coord, ve);
                    Particle pi(coord, vi);

                    // We don't need to use critical section here because each
                    // thread work with him cell
                    electrons->add_particle(pe);
                    ions->add_particle(pi);

                    charged_particle.velocity = v1;

                    // УДАЛЯЕМ нейтрала через swap-and-pop
                    std::swap(neutral_particle, neutrals_data[current_neutral_count - 1]);
                    current_neutral_count--;

                    std::cout << "collision with neutral" << std::endl;
                } else {
                    charged_particle.velocity = v1;
                    neutral_particle.velocity = v2;
                }
            }

            // ФИНАЛЬНОЕ УДАЛЕНИЕ: обрезаем вектор до актуального размера
            if (current_neutral_count < static_cast<int>(neutrals_data.size())) {
                neutrals_data.resize(current_neutral_count);
            }
        }
        if (omp_get_thread_num() == 0) {
            colliderWithNeutrals.profiler.print_report(std::cout);
            colliderWithNeutrals.profiler.reset();
        }
    }
}
