#include "collision.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "Vec.h"

void BinaryCollider::collide_particles(std::vector<ParticlesArray> &species,
                                       int NumPartPerCell, double dt ) {
    for (auto &sp : species) {
        #pragma omp parallel
        for (int j = 0; j < sp.size(); ++j) {
            //std::cout << sp.particlesData(j).size(), sp.charge, density, sp.mass(),
            //    dt << "\n"; 
                double density =
                    sp.particlesData(j).size() / (double) NumPartPerCell;
            collide_same_sort(sp.particlesData(j), sp.charge, density,
                              sp.mass(), dt);
        }
    }

    auto pairs = generatePairs(species.size());
    for (auto &pair : pairs) {
        auto &species1 = species[pair.first];
        auto &species2 = species[pair.second];
#pragma omp parallel
        for (int j = 0; j < species1.size(); ++j) {
            double density1 =
                species1.particlesData(j).size() / (double) NumPartPerCell;
            double density2 =
                species2.particlesData(j).size() / (double) NumPartPerCell;
            collide_diff_sort(species1.particlesData(j), species1.charge,
                              density1, species1.mass(),
                              species2.particlesData(j), species2.charge,
                              density2, species2.mass(), dt);
        }
    }
}

void BinaryCollider::collide_two_particles(double3 &v1, double3 &v2, double q1, double q2,
                                 double n1, double n2, double m1, double m2,
                                 double dt, double variance_factor) {
    const double n = std::min(n1, n2);
    const double m = get_center_mass(m1, m2);
    const double3 u = v1 - v2;
    const double modu = mag(u);
    // double vold1 = v1.x();
    const double variance =
        variance_factor * get_variance_coll(modu, q1, q2, n, m, dt);

    const double sigma = (variance < 1) ? ranomdGenerator.Gauss(sqrt(variance))
                                        : M_PI * ranomdGenerator.Uniform01();

    const double phi = 2 * M_PI * ranomdGenerator.Uniform01();
    const double cosp = cos(phi);
    const double sinp = sin(phi);
    const double sint = 2 * sigma / (1 + sigma * sigma);
    const double cost = 1 - 2 * sigma * sigma / (1 + sigma * sigma);
    const double up = sqrt(u.x() * u.x() + u.y() * u.y());
    double dux, duy, duz;
    if (up < 1.e-8) {
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

void BinaryCollider::collide_same_sort(std::vector<Particle> &particles,
                                              double q1, double n1, double m1,
                                              const double dt) {
    if (particles.size() < 2) {
        return;
    }
    std::vector<int> numbers(particles.size());
    get_numbers_of_colliding_particles(numbers);
    const int startParticle = (numbers.size() % 2 == 0) ? 0 : 3;

    for (size_t i = startParticle; i < numbers.size(); i += 2) {
        double3 v1 = particles[numbers[i]].velocity;
        double3 v2 = particles[numbers[i + 1]].velocity;
        const double variance_factor = 1.0;
        collide_two_particles(v1, v2, q1, q1, n1, n1, m1, m1, dt,
                              variance_factor);
        particles[numbers[i]].velocity = v1;
        particles[numbers[i + 1]].velocity = v2;
    }
    if (numbers.size() % 2 == 1) {
        double3 v1 = particles[numbers[0]].velocity;
        double3 v2 = particles[numbers[1]].velocity;
        double3 v3 = particles[numbers[2]].velocity;
        const double variance_factor = 0.5;
        collide_two_particles(v1, v2, q1, q1, n1, n1, m1, m1, dt,
                              variance_factor);
        collide_two_particles(v1, v3, q1, q1, n1, n1, m1, m1, dt,
                              variance_factor);
        collide_two_particles(v2, v3, q1, q1, n1, n1, m1, m1, dt,
                              variance_factor);
        particles[numbers[0]].velocity = v1;
        particles[numbers[1]].velocity = v2;
        particles[numbers[2]].velocity = v3;
    }
}

void BinaryCollider::collide_diff_sort(std::vector<Particle> &particles1, double q1,
                              double n1, double m1,
                              std::vector<Particle> &particles2, double q2,
                              double n2, double m2, const double dt) {
    const double variance_factor = 1.;
    const int maxSize = std::max(particles1.size(), particles2.size());
    const int minSize = std::min(particles1.size(), particles2.size());

    if (minSize < 1) {
        return;
    }

    std::vector<int> numbers(maxSize);
    get_numbers_of_colliding_particles(numbers);
    for (size_t i = 0; i < numbers.size(); i++) {
        int i1 = numbers[i];
        int i2 = i % minSize;
        if (particles1.size() < particles2.size()){
            std::swap(i1,i2);
        }
        double3 v1 = particles1[i1].velocity;
        double3 v2 = particles2[i2].velocity;
        collide_two_particles(v1, v2, q1, q2, n1, n2, m1, m2, dt,
                              variance_factor);
        particles1[i1].velocity = v1;
        particles2[i2].velocity = v2;
    }
}
