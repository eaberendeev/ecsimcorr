#pragma once

#include <string>

#include "grid/domain.h"
#include "particle_array.h"

class CollisionOperator {
   public:
    virtual ~CollisionOperator() = default;
    virtual void apply(Species& species, const Domain& domain, double dt) = 0;
    virtual std::string info() const = 0;
};
