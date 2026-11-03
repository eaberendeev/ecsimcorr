#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "grid/domain.h"
#include "particle_array.h"

class CollisionOperator {
   public:
    virtual ~CollisionOperator() = default;
    virtual void apply(Species& species, const Domain& domain, double dt) = 0;
    virtual std::string info() const = 0;
};

// Owns the selected processes in their physical application order, including a
// legacy process when chosen by configuration. RNG state survives between steps.
class CollisionManager {
   public:
    CollisionManager() = default;
    explicit CollisionManager(std::vector<std::unique_ptr<CollisionOperator>> processes)
        : processes_(std::move(processes)) {
    }

    void apply(Species& species, const Domain& domain, double dt) {
        for (auto& process : processes_) process->apply(species, domain, dt);
    }
    bool enabled() const {
        return !processes_.empty();
    }

   private:
    std::vector<std::unique_ptr<CollisionOperator>> processes_;
};
