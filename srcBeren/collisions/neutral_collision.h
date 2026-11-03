#pragma once

#include "collisions/collision_manager.h"
#include "random.h"

struct NeutralCollisionConfig {
    std::string charged_species;
    std::string neutral_species;
    std::string electron_product;
    std::string ion_product;
    std::string scheme = "physical_only";
    bool electron_ionization = true;
    bool proton_ionization = true;
    bool proton_charge_exchange = true;
};

// Collisions change velocities, populations and products within the current cell.
// Particle transport and boundary conditions belong to the simulation step.
class NeutralCollisionOperator : public CollisionOperator {
   public:
    NeutralCollisionOperator(const NeutralCollisionConfig& config, double n0);
    void apply(Species& species, double dt) override;
    std::string info() const override;

   private:
    NeutralCollisionConfig config_;
    double n0_;
    LehmerEngine baseRndEng_;
};
