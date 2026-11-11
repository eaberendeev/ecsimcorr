#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#include "drift_kinetic_push.h"
#include "particles/Particle.h"
#include "utils/containers.h"

// Change of the model of one marker between the kinetic description
// (position x, velocity v) and the guiding-centre one (R, v_par, mu).
//
// The marker keeps its position, R = x. The charge density of both models uses
// the same shape, so the change does not touch the charge on the grid and
// Gauss's law is kept exactly. The price is a displacement of the centre by
// the Larmor radius rho, which the criterion keeps below a fraction of a cell.
//
// The energy is kept exactly: m v^2 / 2 = m v_par^2 / 2 + mu |B(R)|, with the
// same |B| that the guiding-centre energy diagnostic uses. The whole
// perpendicular energy goes into mu, including the E x B drift energy, which
// the guiding-centre model does not carry separately. The perpendicular
// momentum of the gyration is not transferred: the guiding centre has none.
//
// Which model fits a marker. The kinetic orbit of the implicit midpoint rule
// turns by the angle theta per step, tan(theta / 2) = omega_c dt / 2, on a
// circle of the exact Larmor radius, but the current it deposits has the
// magnetic moment mu / (1 + (omega_c dt / 2)^2) instead of mu. For
// omega_c dt >> 1 a kinetic marker carries almost no diamagnetic current,
// while a guiding centre carries all of it as the magnetization -mu b. The two
// descriptions agree only while omega_c dt is small, so the model is chosen by
// omega_c dt = |q| |B| dt / m at the marker: above an upper threshold a
// guiding centre, below a lower one kinetic (`numerical_moment_factor`). This
// rule depends on the position only: it does not split the velocity
// distribution at one place into two parts with different responses.
//
// The two models also answer differently to fields at frequencies near
// omega_c, which the kinetic model resolves when omega_c dt < 1: a kinetic
// marker has inertia across B, a guiding centre follows E x B at once. A
// change of model therefore replaces the response of the marker to such
// fields, which are mostly particle noise, and adds to that noise. Markers
// that cross a border between regions one by one add little; a whole region
// changing its model at once does not (measured in docs/drift-kinetic.md).
// In addition a guiding centre must have its Larmor radius below a given
// number of cells: this bounds the displacement R = x and the gradient of B
// across the orbit.
//
// The decision must not depend on the gyrophase. In an electric field the
// speed across B oscillates over a gyration between |u - v_E| and |u + v_E|
// (u: gyration speed, v_E = E x B / B^2). A criterion on the instantaneous
// speed would pick the markers at the slow phase, i.e. on one side of their
// Larmor circles: with R = x the centres would shift along E, each change of
// model would polarize the plasma and feed the field. The Larmor radius of a
// kinetic marker is therefore taken in the frame of the E x B drift, where it
// is constant over the gyration, like mu of a guiding centre.
//
// All quantities belong to the same integer time level n.
namespace driftKinetic {

// Larmor radius m |v_perp - v_E| / (|q| |B|) of a kinetic marker in the fields
// E, B at its position; infinite without a magnetic field.
inline double larmor_radius(const Vector3R& velocity, const Vector3R& electric, const Vector3R& magnetic, double charge,
                            double mass) {
    const double lenB = magnetic.length();
    if (!(lenB > 0.0) || charge == 0.0)
        return std::numeric_limits<double>::infinity();
    const Vector3R b = magnetic / lenB;
    const Vector3R gyration = velocity - velocity.dot(b) * b - electric.cross(b) / lenB;
    return mass * gyration.length() / (std::abs(charge) * lenB);
}

// Ratio of the magnetic moment of the current deposited by a kinetic orbit of
// the implicit midpoint rule to the moment m v_perp^2 / (2 |B|) of its energy.
inline double numerical_moment_factor(double omega_dt) {
    return 1.0 / (1.0 + 0.25 * omega_dt * omega_dt);
}

// Larmor radius of a guiding centre, v_perp = sqrt(2 mu |B| / m).
inline double larmor_radius(const GuidingCenter& centre, const Vector3R& magnetic) {
    const double lenB = magnetic.length();
    if (!(lenB > 0.0) || centre.charge == 0.0)
        return std::numeric_limits<double>::infinity();
    return std::sqrt(2.0 * centre.mu * centre.mass / lenB) / std::abs(centre.charge);
}

// Kinetic marker -> guiding centre in the field `magnetic` = B(x), |B| > 0.
inline GuidingCenter to_guiding_center(const Particle& particle, const Vector3R& magnetic, double charge, double mass,
                                       double weight) {
    const double lenB = magnetic.length();
    const Vector3R b = magnetic / lenB;
    GuidingCenter centre;
    centre.position = particle.coord;
    centre.v_parallel = particle.velocity.dot(b);
    // v_perp^2 = v^2 - v_par^2: the energy identity holds to rounding.
    const double perp2 = std::max(0.0, particle.velocity.squaredNorm() - centre.v_parallel * centre.v_parallel);
    centre.mu = 0.5 * mass * perp2 / lenB;
    centre.charge = charge;
    centre.mass = mass;
    centre.weight = weight;
    return centre;
}

// Guiding centre -> kinetic marker in the field `magnetic` = B(R), |B| > 0.
// `gyrophase` (radians) fixes the direction of v_perp in the plane across b;
// the guiding centre does not know it. The marker gets no E x B drift
// velocity: its whole perpendicular energy is gyration, as mu holds it. In a
// slowly varying E this shifts the centre of its Larmor circle by
// v_E x b / omega_c, the polarization displacement of a particle released at
// rest. The drift is not added because the field at one time level is mostly
// particle noise at frequencies above omega_c, for which E x B / B^2 is not
// the velocity of the plasma.
inline Particle to_kinetic(const GuidingCenter& centre, const Vector3R& magnetic, double gyrophase) {
    const double lenB = magnetic.length();
    const Vector3R b = magnetic / lenB;
    // Orthonormal pair across b, built from the axis least aligned with b.
    Vector3R axis(1, 0, 0);
    if (std::abs(b.y()) <= std::abs(b.x()) && std::abs(b.y()) <= std::abs(b.z()))
        axis = Vector3R(0, 1, 0);
    else if (std::abs(b.z()) <= std::abs(b.x()) && std::abs(b.z()) <= std::abs(b.y()))
        axis = Vector3R(0, 0, 1);
    const Vector3R e1 = b.cross(axis).normalized();
    const Vector3R e2 = b.cross(e1);
    const double v_perp = std::sqrt(2.0 * centre.mu * lenB / centre.mass);
    Particle particle;
    particle.coord = centre.position;
    particle.velocity = centre.v_parallel * b + v_perp * (std::cos(gyrophase) * e1 + std::sin(gyrophase) * e2);
    return particle;
}

// Gyrophase in [0, 2 pi) as a hash of the marker state and the time step.
// It does not depend on the thread count or on the order of the markers, and
// a restarted run repeats it.
inline double hashed_gyrophase(const GuidingCenter& centre, int timestep) {
    const auto mix = [](std::uint64_t z) {   // splitmix64 finalizer
        z += 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    };
    std::uint64_t state = mix(static_cast<std::uint64_t>(timestep));
    for (const double value : {centre.position.x(), centre.position.y(), centre.position.z(), centre.v_parallel}) {
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        state = mix(state ^ bits);
    }
    return 2.0 * M_PI * static_cast<double>(state >> 11) / 9007199254740992.0;   // 53 random bits
}

}   // namespace driftKinetic
