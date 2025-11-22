//#include "Particles.h"

#include "collision.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "Vec.h"
#include "collisions_with_neutrals.h"

double BinaryCollider::get_variance_coll(double u, double q1, double q2,
                                                double n, double m, double dt) {
    const double lk = 15;
    return (pow(SGS::get_plasma_freq(n0), 3) /
            (SGS::c * SGS::c * SGS::c * n0)) *
           (lk * q1 * q1 * q2 * q2 * n * dt) / (8 * M_PI * m * m * u * u * u);
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
  const int firstGroup =
      r * (q + 1); // number of particles in first group in sort1

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

void BinaryCollider::bin_collide(double3 &v1, double3 &v2, double q1, double q2, double n1,
                 double n2, double m1, double m2, double dt,
                 double variance_factor) {
  const double n = std::min(n1, n2);
  const double m = get_center_mass(m1, m2);
  const double3 u = v1 - v2;
  const double modu = mag(u);
  const double variance =
      variance_factor * get_variance_coll(modu, q1, q2, n, m, dt);

  const double sigma =
      (variance < 1) ? gen.Gauss(sqrt(variance)) : M_PI * gen.Uniform01();

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
    dux = (u.x() / up) * u.z() * sint * cosp -
          (u.y() / up) * modu * sint * sinp - u.x() * (1 - cost);
    duy = (u.y() / up) * u.z() * sint * cosp +
          (u.x() / up) * modu * sint * sinp - u.y() * (1 - cost);
    duz = -up * sint * cosp - u.z() * (1 - cost);
  }
  const double3 du = double3(dux, duy, duz);
  v1 += (m / m1) * du;
  v2 -= (m / m2) * du;
}

void BinaryCollider::collide_same_sort_binary(Species &species, const double dt) {
    for (auto &sp : species) {
        const double q = sp->charge;
        const double m1 = sp->mass();
#pragma omp parallel for schedule(dynamic, 32)
        for (auto pk = 0; pk < sp->size(); pk++) {
            BinaryCollisionSameType collider(sp->particlesData(pk).size(),
                                             gen.gen());
            while (collider.canCollide()) {
                auto pair = collider.get_pair();
                double3 v1 = sp->particlesData(pk)[pair.first].velocity;
                double3 v2 = sp->particlesData(pk)[pair.second].velocity;
                double n1 =
                    sp->particlesData(pk).size() / (double) sp->NumPartPerCell;
                double variance_factor = collider.get_variance_factor();
                bin_collide(v1, v2, q, q, n1, n1, m1, m1, dt, variance_factor);
                sp->particlesData(pk)[pair.first].velocity = v1;
                sp->particlesData(pk)[pair.second].velocity = v2;
            }
        }
    }
}

void BinaryCollider::collide_ion_electron_binary(
    Species &species, const double dt) {
    int electrons = get_num_of_type_particles(species, "Electrons");
    int ions = get_num_of_type_particles(species, "Ions");
    const double q1 = species[electrons]->charge;
    //const int n1 = species[electrons].density;
    const double m1 = species[electrons]->mass();
    const double q2 = species[ions]->charge;
    //const int n2 = species[ions].density;
    const double m2 = species[ions]->mass();
#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < species[electrons]->size(); pk++) {
        BinaryCollisionDiffType collider(
            species[electrons]->particlesData(pk).size(),
            species[ions]->particlesData(pk).size(), gen.gen());
        while (collider.canCollide()) {
            auto pair = collider.get_pair();
            double3 v1 =
                species[electrons]->particlesData(pk)[pair.first].velocity;
            double3 v2 = species[ions]->particlesData(pk)[pair.second].velocity;
            const double variance_factor = 1.;
            double n1 = species[electrons]->particlesData(pk).size() /
                        (double) species[electrons]->NumPartPerCell;
            double n2 = species[ions]->particlesData(pk).size() /
                        (double) species[ions]->NumPartPerCell;
            bin_collide(v1, v2, q1, q2, n1, n2, m1, m2, dt, variance_factor);
            species[electrons]->particlesData(pk)[pair.first].velocity = v1;
            species[ions]->particlesData(pk)[pair.second].velocity = v2;
        }
    }
}

void BinaryCollider::collide_with_neutrals_binary(Species &species,
                                                  const double dt) {
    for (auto& sp : species) {
        sp->update_count_in_cell();
    }

    for (auto i = 0; i < species.size(); i++) {
        if (species[i]->is_neutral() )
            continue;
        collide_with_neutrals_binary_impl(species, i, dt);
        species[i]->update_count_in_cell();
    }
}

void BinaryCollider::collide_with_neutrals_binary_impl(Species &species,
                                                       const int pType,
                                                       const double dt) {
    int neutrals_type = get_num_of_type_particles(species, "Neutrals");
    int electrons = get_num_of_type_particles(species, "Electrons");
    int ions = get_num_of_type_particles(species, "Ions");
    if (neutrals_type < 0 || electrons < 0 || ions < 0) {
      std::cout << "Error: Invalid type of species for collision\n";
      return;
    }
    const double m1 = species[pType]->mass();
    const double m2 = species[neutrals_type]->mass();

#pragma omp parallel 
{
    ColliderWithNeutrals colliderWithNeutrals(n0);

#pragma omp for schedule(dynamic, 32)
    for (auto pk = 0; pk < species[pType]->size(); pk++) {

        int pInCell = species[pType]->countInCell(pk);
        int nInCell = species[neutrals_type]->countInCell(pk);

        if (pInCell == 0 || nInCell == 0)
            continue;

        double n1 = pInCell / (double) species[pType]->NumPartPerCell;
        double n2 = nInCell / (double) species[neutrals_type]->NumPartPerCell;

        // Работаем напрямую с данными нейтралов
        auto &neutrals_data = species[neutrals_type]->particlesData(pk);
        int current_neutral_count = nInCell;

        for (int i = 0; i < pInCell && current_neutral_count > 0; i++) {
            std::uniform_int_distribution<> dis(0, current_neutral_count - 1);
            int randomIndex = dis(gen.gen());

            Particle &charged_particle = species[pType]->particlesData(pk)[i];
            Particle &neutral_particle = neutrals_data[randomIndex];

            double3 v1 = charged_particle.velocity;
            double3 v2 = neutral_particle.velocity;

            auto [is_collided, ve, vi] =
                colliderWithNeutrals.collision_with_neutral(v1, v2, m1, m2, n1,
                                                            n2, dt, 1. / dt);
            if (is_collided) {
                double3 coord = neutral_particle.coord;
                Particle pe(coord, ve);
                Particle pi(coord, vi);

                // We don't need to use critical section here because each thread work with 
                // him cell
                species[electrons]->add_particle(pe);
                species[ions]->add_particle(pi);
                
                charged_particle.velocity = v1;

                // УДАЛЯЕМ нейтрала через swap-and-pop (как в вашем коде)
                std::swap(neutral_particle,
                          neutrals_data[current_neutral_count - 1]);
                current_neutral_count--;

                std::cout << "collision with neutral" << std::endl;
            } else {
                charged_particle.velocity = v1;
                neutral_particle.velocity = v2;
            }
        }

        // ФИНАЛЬНОЕ УДАЛЕНИЕ: обрезаем вектор до актуального размера
        if (current_neutral_count < neutrals_data.size()) {
            neutrals_data.resize(current_neutral_count);
        }
    }
    if (omp_get_thread_num() == 0) {
        colliderWithNeutrals.profiler.print_report(std::cout);
        colliderWithNeutrals.profiler.reset();
    }
  }

    species[neutrals_type]->update_count_in_cell();
}