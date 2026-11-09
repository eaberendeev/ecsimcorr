#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "drift_kinetic_push.h"

namespace driftKinetic {

// Recursively halves the guiding-centre step when the implicit midpoint
// Picard/Aitken solve does not converge.  A failed full step is retried as two
// half steps; the accepted trajectory is the concatenation of the leaf steps.
class AdaptiveSubstepping {
   public:
    void configure(bool enabled, int max_depth) {
        if (max_depth < 0 || max_depth > 30)
            throw std::invalid_argument("drift-kinetic max_substep_depth must be in [0, 30]");
        enabled_ = enabled;
        max_depth_ = max_depth;
    }

    int max_depth() const {
        return max_depth_;
    }

    struct Result {
        GuidingCenter state;
        int leaves = 0;
        int depth = 0;
        bool converged = false;
    };

    Result push(const GuidingCenterPush& pusher, const GuidingCenter& initial, double dt,
                const FieldFunction& fields) const {
        Result result;
        if (!enabled_ || max_depth_ == 0) {
            const auto one = pusher.push(initial, dt, fields, true);
            result.state = one.state;
            result.leaves = 1;
            result.converged = one.converged;
            return result;
        }
        result = push_recursive(pusher, initial, dt, fields, 0);
        return result;
    }

   private:
    Result push_recursive(const GuidingCenterPush& pusher, const GuidingCenter& initial, double dt,
                          const FieldFunction& fields, int depth) const {
        const auto one = pusher.push(initial, dt, fields, true);
        if (one.converged || depth >= max_depth_) {
            Result result;
            result.state = one.state;
            result.leaves = 1;
            result.depth = depth;
            result.converged = one.converged;
            return result;
        }

        const auto left = push_recursive(pusher, initial, 0.5 * dt, fields, depth + 1);
        const auto right = push_recursive(pusher, left.state, 0.5 * dt, fields, depth + 1);
        Result result;
        result.state = right.state;
        result.leaves = left.leaves + right.leaves;
        result.depth = std::max(left.depth, right.depth);
        result.converged = left.converged && right.converged;
        return result;
    }

    bool enabled_ = false;
    int max_depth_ = 4;
};

}  // namespace driftKinetic
