#include "collisions/collision_manager.h"

void CollisionManager::apply(Species& species, const Domain& domain, double dt) {
    if (operators_.empty()) {
        if (legacy_fallback_)
            legacy_fallback_->apply(species, domain, dt);
        return;
    }
    for (auto& op : operators_) {
        op->apply(species, domain, dt);
    }
}
