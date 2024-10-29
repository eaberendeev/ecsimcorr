//#include "Particles.h"
#ifndef COLLISION_OLD_H_
#define COLLISION_OLD_H_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>
#include "sgs.h"
#include "ParticlesArray.h"

class BinaryCollider {
  public:
    double get_variance_coll(double u, double q1, double q2, double n,
                                   double m, double dt);
       BinaryCollider(double n0)
       : n0(n0) {
        gen.SetRandSeed(13);
   }
   double n0;

  public:
   void collide_same_sort_binary(Species &species,
                                 const double dt);
   void collide_ion_electron_binary(Species &species, const double dt);
   void bin_collide(double3 &v1, double3 &v2, double q1, double q2, double n1,
                    double n2, double m1, double m2, double dt,
                    double variance_factor);

   ThreadRandomGenerator gen;
};
// Takizuka. A Binary Collision Model for Plasma Simulation
// with a Particle Code // JOURNAL OF COMPUTATIONAL PHYSICS 25, 205-219
// (1977)

// Takizuka1977 / case 1. same type of particles
class BinaryCollisionSameType {
   public:
    BinaryCollisionSameType(int size, std::mt19937 &g) : v(size) {
        //std::mt19937 g;
        for (size_t i = 0; i < v.size(); i++) {
            v[i] = i;
        }
        std::shuffle(v.begin(), v.end(), g);
        // std::copy(v.begin(), v.end(),
        //           std::ostream_iterator<int>(std::cout, " "));
        // std::cout << '\n';
        event = -1;
        isOddEvents = (size % 2 == 1) ? true : false;
        maxEvent = size / 2 + 2 * (size % 2);
    }
    bool canCollide() { return (v.size() > 1) && (event < maxEvent - 1); }
    double get_variance_factor() {
        // Takizuka1977 / case 1b
        if (!isOddEvents || event >= 3) {
            return 1.;
        }
        return 0.5;
    }

    std::pair<int, int> get_pair();

   private:
    std::vector<int> v;
    int event;
    int maxEvent;
    bool isOddEvents;   // check for case 1a or 1b
};

// Takizuka1977 / case 2. different type of particles
class BinaryCollisionDiffType {
public:
  BinaryCollisionDiffType(int size1, int size2, std::mt19937 &g)
      : v1(std::max(size1, size2)), v2(std::min(size1, size2)) {
    for (size_t i = 0; i < v1.size(); i++) {
      v1[i] = i;
    }
    for (size_t i = 0; i < v2.size(); i++) {
      v2[i] = i;
    }
    std::shuffle(v1.begin(), v1.end(), g);
    std::shuffle(v2.begin(), v2.end(), g);
    // std::copy(v1.begin(), v1.end(), std::ostream_iterator<int>(std::cout, " "));
    // std::cout << '\n';
    // std::copy(v2.begin(), v2.end(), std::ostream_iterator<int>(std::cout, " "));
    // std::cout << '\n';
    event = -1;
    maxEvent = std::max(size1, size2);
    leadingSort = (size1 >= size2) ? 0 : 1;
    firstInd = 0;
  }
  bool canCollide() {
    return (v1.size() == v2.size() || (v1.size() > 1 && v2.size() > 0)) &&
           (event < maxEvent - 1);
  }
  std::pair<int, int> get_pair();

private:
  std::vector<int> v1, v2;
  int leadingSort; // sort with large number of particles
  int event;
  int maxEvent;
  int firstInd; // current index of first array
};

inline double get_center_mass(double m1, double m2) {
  return m1 * m2 / (m1 + m2);
}


#endif   // COLLISION_H_