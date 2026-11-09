#pragma once

#include "grid/domain.h"
#include "particles/Particle.h"
#include "utils/containers.h"

// Cubical Whitney forms: nodal charge, edge E, face B; normalized units c = 1.
namespace symplectic {
void fill_periodic_ghosts(const Domain& domain, Field3d& field);

template <bool isMagnetic>
Vector3R gather(const Domain& domain, const Field3d& field, const Vector3R& position);

void deposit_charge(const Domain& domain, Field3d& rho, const Vector3R& position, double weight);

void electric_flow(const Domain& domain, Field3d& B, const Field3d& E, double dt);

void magnetic_flow(const Domain& domain, Field3d& E, const Field3d& B, double dt);

// Requires abs(dt * v_axis) <= cell_size(axis); processes at most two cells.
void advance_axis(const Domain& domain, Particle& p, int axis, double dt, double qm, double qw, const Field3d& B,
                  Field3d& integral);

void drift_neutral(const Domain& domain, Particle& particle, double dt);

double divergence(const Domain& domain, const Field3d& field, int i, int j, int k, bool isMagnetic);

double field_energy(const Domain& domain, const Field3d& field);

// Execute Hx,Hy,Hz (or the reverse order) per particle, sharing a bucket-local field/current patch.
void advance_bucket(const Domain& domain, std::vector<Particle>& particles, int bucket, bool reverse, double dt,
                    double qm, double qw, const Field3d& magnetic, Field3d& integral);
}   // namespace symplectic
