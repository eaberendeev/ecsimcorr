#ifndef PARTICLE_H_
#define PARTICLE_H_
#include <assert.h>

#include <functional>

#include "containers.h"
#include "sgs.h"

struct ParticleBase {
    Vector3R coord;
    Vector3R velocity;

    ParticleBase() = default;
    ParticleBase(double x, double y, double z, double vx, double vy, double vz) : coord{x, y, z}, velocity{vx, vy, vz} {
    }
    ParticleBase(const Vector3R& x, const Vector3R& v) : coord(x), velocity(v) {
    }

    void move(double dt) {
        coord += velocity * dt;
    }
};

struct ParticleSimple : public ParticleBase {
    Vector3R initCoord;
    Vector3R initVelocity;

    ParticleSimple() = default;

    ParticleSimple(double x, double y, double z, double vx, double vy, double vz)
        : ParticleBase(x, y, z, vx, vy, vz), initCoord{x, y, z}, initVelocity{vx, vy, vz} {
    }

    ParticleSimple(const Vector3R& x, const Vector3R& v) : ParticleBase(x, v), initCoord(x), initVelocity(v) {
    }

#ifdef SET_PARTICLE_IDS
    size_t id;
#endif

    friend std::ostream& operator<<(std::ostream& out, const ParticleSimple& particle);
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
    friend std::ostream& operator<<(std::ostream& out, const ParticleMPW& particle);
};
struct ParticleMass : ParticleSimple {
    double mass;
    friend std::ostream& operator<<(std::ostream& out, const ParticleMass& particle);
};

struct ParticleMassMPW : ParticleSimple {
    double mass, mpw;
    friend std::ostream& operator<<(std::ostream& out, const ParticleMassMPW& particle);
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
inline double get_energy_particle(const Vector3R& v, const double m, const double mpw) {
    return 0.5 * mpw * m * v.dot(v);
}

inline double get_energy_particle(const double v, const double m, const double mpw) {
    return 0.5 * mpw * m * v * v;
}

inline double get_energy_particle(const double v1, const double v2, const double m, const double mpw) {
    return 0.5 * mpw * m * (v1 * v1 + v2 * v2);
}

inline double get_energy_particles(const std::vector<Particle>& particles, double mass, double mpw) {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (const auto& particle : particles) {
        energy += get_energy_particle(particle.velocity, mass, mpw);
    }

    return energy;
}
inline double convert_kev_to_velocity(double kev, double mass) {
    return sqrt(2 * kev / SGS::MC2 / mass);
}

inline Vector3R convert_kev_to_sigma(Vector3R sigma, double mass) {
    double sigmax = std::sqrt(sigma.x() / SGS::MC2 / mass);
    double sigmay = std::sqrt(sigma.y() / SGS::MC2 / mass);
    double sigmaz = std::sqrt(sigma.z() / SGS::MC2 / mass);
    return Vector3R(sigmax, sigmay, sigmaz);
}

#endif
