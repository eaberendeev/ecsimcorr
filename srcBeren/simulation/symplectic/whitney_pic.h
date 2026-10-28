#pragma once

#include "Particle.h"
#include "World.h"

// Lowest-order cubical Whitney complex: nodal charge, edge E, face B.
// All field coefficients use Beren's normalized units (c = 1). Only the
// num_cells unique periodic nodes are evolved; ghosts are copies for outputs.
namespace symplectic {
class WhitneyPIC {
   public:
    explicit WhitneyPIC(const Domain& domain) : domain_(domain) {
    }
    void sync(Field3d& field) const;
    Vector3R gather(const Field3d& field, const Vector3R& position, bool magnetic) const;
    void deposit_charge(Field3d& rho, const Vector3R& position, double charge_density_weight) const;
    void electric_flow(Field3d& B, const Field3d& E, double dt) const;
    void magnetic_flow(Field3d& E, const Field3d& B, double dt) const;
    // Exact H_axis flow for one particle. current_integral accumulates integral J dt
    // over the COMPLETE step, including both half flows. Sequential shared deposits.
    void kinetic_flow(Particle& particle, int axis, double dt, double charge_mass, double charge_density_weight,
                      Field3d& E, const Field3d& B, Field3d& current_integral) const;
    double divergence(const Field3d& field, int i, int j, int k, bool magnetic) const;
    double field_energy(const Field3d& field) const;

   private:
    const Domain& domain_;
    int index(int i, int axis) const;
    double value(const Field3d& field, Vector3I node, int component) const;
};
}   // namespace symplectic
