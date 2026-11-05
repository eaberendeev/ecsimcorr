#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "particles/particle_array.h"

class CollisionOperator {
   public:
    virtual ~CollisionOperator() = default;
    virtual void apply(Species& species, double dt) = 0;
    virtual std::string info() const = 0;
};

// Owns the configured processes in physical application order.
// RNG state survives between steps; processes do not transport particles.
class CollisionManager {
   public:
    CollisionManager() = default;
    explicit CollisionManager(std::vector<std::unique_ptr<CollisionOperator>> processes)
        : processes_(std::move(processes)) {
    }

    void apply(Species& species, double dt) {
        for (auto& process : processes_) process->apply(species, dt);
    }
    bool enabled() const {
        return !processes_.empty();
    }

   private:
    std::vector<std::unique_ptr<CollisionOperator>> processes_;
};
