#ifndef PARTICLE_H_
#define PARTICLE_H_
#include <assert.h>

#include <functional>

#include "containers.h"

#define SET_PARTICLE_IDS

struct ParticleSimple {
    double3 coord;
    double3 velocity;
    double3 initCoord;
    double3 initVelocity;
    ParticleSimple() {}
    ParticleSimple(double x,double y, double z, double vx, double vy, double vz){
        coord = {x, y, z};
        velocity = {vx, vy, vz};
        initCoord = {x, y, z};
        initVelocity = {vx, vy, vz};
    }
    ParticleSimple(const double3& x, const double3& v) {
        coord = x;
        velocity = v;
        initCoord = x;
        initVelocity = v;
    }
#ifdef SET_PARTICLE_IDS
    size_t id;
#endif

    friend std::ostream& operator<<(std::ostream& out,
                                    const ParticleSimple& particle);

    void move(double dt) { coord += velocity * dt; }
};

inline int pos_ind(int index, int n, int _size1, int _size2, int _size3) {
    std::vector<int> dim = {_size1, _size2, _size3};
    int capacity = 1;
    for (unsigned int i = n + 1; i < dim.size(); i++) {
        capacity *= dim[i];
    }
    return (index / capacity) % dim[n];
}

struct ParticleMPW : ParticleSimple {
    double mpw;
    friend std::ostream& operator<<(std::ostream& out,
                                    const ParticleMPW& particle);
};
struct ParticleMass : ParticleSimple {
    double mass;
    friend std::ostream& operator<<(std::ostream& out,
                                    const ParticleMass& particle);
};

struct ParticleMassMPW : ParticleSimple {
    double mass, mpw;
    friend std::ostream& operator<<(std::ostream& out,
                                    const ParticleMassMPW& particle);
};

#ifdef PARTICLE_MASS
#ifdef PARTICLE_MPW
typedef ParticleMassMPW Particle;
#else
typedef ParticleMass Particle;
#endif
#else
#ifdef PARTICLE_MPW
typedef ParticleMPW Particle;
#else
typedef ParticleSimple Particle;
#endif
#endif

/**
 * Calculates the kinetic energy for a particle with the given velocity,
 * mass, and mass per weight. Overloaded for vector, scalar velocity.
 */
inline double get_energy_particle(const double3& v, const double m,
                                  const double mpw) {
    return 0.5 * mpw * m * dot(v, v);
}

inline double get_energy_particle(const double v, const double m,
                                  const double mpw) {
    return 0.5 * mpw * m * v * v;
}

inline double get_energy_particle(const double v1, const double v2,
                                  const double m, const double mpw) {
    return 0.5 * mpw * m * (v1 * v1 + v2 * v2);
}

inline double get_energy_particles(const std::vector<Particle>& particles,
                                  double mass, double mpw) {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (const auto& particle : particles) {
        energy += get_energy_particle(particle.velocity, mass, mpw);
    }

    return energy;
}

#endif
