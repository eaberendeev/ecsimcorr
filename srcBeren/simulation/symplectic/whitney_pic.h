#pragma once

#include "particles/particle_array.h"
#include "utils/containers.h"

// Lowest-order cubical Whitney forms on the Yee mesh: nodal charge, edge E, face B; c = 1.
// The kernels index storage directly, without periodic wrapping: fields must hold
// their periodic images in the ghost nodes, particles must lie inside the box, and
// deposits into ghost nodes are left for the periodic boundary condition to fold.
namespace symplectic {
// v += impulse * E(x) for every particle.
void kick_velocities(ParticlesArray& particles, const Field3d& E, double impulse);

// Exact flows of H_x, H_y, H_z (H_z, H_y, H_x when reverse) of every particle at fixed B.
// H_a moves x_a by dt*v_a, rotates the transverse velocities by the path integral of B and
// adds q*mpw times the path integral of the a-edge forms to `current`. Coordinates are
// not wrapped and buckets are not updated. Every subflow must end in its starting cell
// or in a neighbor, so abs(dt*v_a) may not exceed cell_size(a); otherwise the particle
// is left between subflows and false is returned.
bool kinetic_flow(ParticlesArray& particles, const Field3d& B, double dt, bool reverse, Field3d& current);

// Maximum of abs(div B) over the cells.
double max_magnetic_divergence(const Domain& domain, const Field3d& B);
}   // namespace symplectic
