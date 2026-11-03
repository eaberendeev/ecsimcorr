#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "collisions/collision_operator.h"

// Owns the ordered channels for one simulation. Process RNG state survives steps.
class CollisionManager {
   public:
    void add(std::unique_ptr<CollisionOperator> process) {
        operators_.push_back(std::move(process));
    }
    void set_legacy_fallback(std::unique_ptr<CollisionOperator> process) {
        legacy_fallback_ = std::move(process);
    }
    void apply(Species& species, const Domain& domain, double dt);
    bool has_collisions() const {
        return !operators_.empty();
    }
    bool enabled() const {
        return has_collisions() || bool(legacy_fallback_);
    }
    std::string legacy_info() const {
        return legacy_fallback_ ? legacy_fallback_->info() : "None";
    }

   private:
    std::vector<std::unique_ptr<CollisionOperator>> operators_;
    std::unique_ptr<CollisionOperator> legacy_fallback_;
};
