#ifndef COLLISION_H_
#define COLLISION_H_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "Particle.h"
#include "random_generator.h"
#include "sgs.h"
#include "Vec.h"

inline std::vector<std::pair<int, int>> generatePairs(int n) {
    std::vector<std::pair<int, int>> pairs;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            pairs.push_back(std::make_pair(i, j));
        }
    }

    return pairs;
}

inline double get_center_mass(double m1, double m2) {
    return m1 * m2 / (m1 + m2);
}

// Takizuka. A Binary Collision Model for Plasma Simulation
// with a Particle Code // JOURNAL OF COMPUTATIONAL PHYSICS 25, 205-219
// (1977)
class BinaryCollider {
   public:
    BinaryCollider(double n0, double lk = 15) : n0(n0), lk(lk){};
    void collide_same_sort(std::vector<Particle> & particles, double q1,
                                  double n1, double m1, const double dt);
    void collide_diff_sort(std::vector<Particle> & particles1, double q1,
                                  double n1, double m1,
                                  std::vector<Particle>& particles2, double q2,
                                  double n2, double m2, const double dt);

   private:
    RandomGenerator ranomdGenerator;
    double n0;
    double lk;   // Culon lagarifm
    void collide_two_particles(double3 & v1, double3 & v2, double q1, double q2,
                               double n1, double n2, double m1, double m2,
                               double dt, double variance_factor);
    double get_variance_coll(double u, double q1, double q2, double n, double m,
                             double dt) {
        return (pow(SGS::get_plasma_freq(n0), 3) /
                (SGS::c * SGS::c * SGS::c * n0)) *
               (lk * q1 * q1 * q2 * q2 * n * dt) /
               (8 * M_PI * m * m * u * u * u);
    }
    void get_numbers_of_colliding_particles(std::vector<int> & numbers) {
        for (size_t i = 0; i < numbers.size(); i++) {
            numbers[i] = i;
        }
        std::shuffle(numbers.begin(), numbers.end(), ranomdGenerator.gen);
    }
};

#endif   // COLLISION_H_