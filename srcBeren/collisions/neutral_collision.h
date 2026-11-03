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

class NeutralCollisionOperator : public CollisionOperator {
   public:
    NeutralCollisionOperator(const NeutralCollisionConfig& config, double n0);
    void apply(Species& species, const Domain& domain, double dt) override;
    std::string info() const override;

   private:
    NeutralCollisionConfig config_;
    double n0_;
    LehmerEngine baseRndEng_;
};
