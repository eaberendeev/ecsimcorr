#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "drift_kinetic_push.h"

namespace driftKinetic {

// Fields of the orbit start -> end over the step fraction [theta0, theta1].
using IntervalFieldFunction =
    std::function<Fields(const Vector3R& start, const Vector3R& end, double theta0, double theta1)>;

// Counters of one field evaluation; merged over particles and threads.
struct PushStatistics {
    long particles = 0;
    long root_failures = 0;      // whole-step solves that did not converge
    long failed_particles = 0;   // particles with a non-converged leaf
    long leaves = 0;
    long iteration_sum = 0;
    long retries = 0;         // interval bisections
    long leaf_failures = 0;   // accepted leaves that did not converge
    int max_iterations = 0;
    int depth = 0;
    double max_leaf_residue_r = 0.0;
    double max_leaf_residue_v = 0.0;

    double average_iterations() const {
        return particles ? static_cast<double>(iteration_sum) / particles : 0.0;
    }

    void merge(const PushStatistics& other) {
        particles += other.particles;
        root_failures += other.root_failures;
        failed_particles += other.failed_particles;
        leaves += other.leaves;
        iteration_sum += other.iteration_sum;
        retries += other.retries;
        leaf_failures += other.leaf_failures;
        max_iterations = std::max(max_iterations, other.max_iterations);
        depth = std::max(depth, other.depth);
        max_leaf_residue_r = std::max(max_leaf_residue_r, other.max_leaf_residue_r);
        max_leaf_residue_v = std::max(max_leaf_residue_v, other.max_leaf_residue_v);
    }
};

// One accepted piece [begin, end] (fractions of dt) of a guiding-centre step.
struct Leaf {
    double begin = 0.0, end = 1.0;
    GuidingCenter start, finish;
};

// Energy-consistent adaptive substepping (xpic
// drift_kinetic::AdaptiveSubstepping). A whole-step solve that does not
// converge is retried as two half intervals, recursively. Every leaf is
// pushed with the magnetic field interpolated to its own time interval and
// deposits its current and magnetization with the weight (end - begin), so
// the sum over leaves is the step source used by Maxwell's equations.
class AdaptiveSubstepping {
   public:
    void configure(bool enabled, int max_depth) {
        if (max_depth < 0 || max_depth > 30)
            throw std::invalid_argument("drift-kinetic max_substep_depth must be in [0, 30]");
        enabled_ = enabled;
        max_depth_ = max_depth;
    }

    int max_depth() const {
        return enabled_ ? max_depth_ : 0;
    }

    // Returns false for a non-finite orbit; `leaves` then must not be deposited.
    // `guess` is the end state of the previous field evaluation.
    bool advance(const GuidingCenterPush& pusher, const GuidingCenter& initial, const GuidingCenter& guess, double dt,
                 const IntervalFieldFunction& fields, bool use_aitken, std::vector<Leaf>& leaves,
                 PushStatistics& stats) const {
        leaves.clear();
        PushStatistics step;
        step.particles = 1;
        GuidingCenter finish = guess;
        const bool valid = advance(pusher, 0.0, 1.0, initial, finish, 0, dt, fields, use_aitken, leaves, step);
        if (valid) {
            step.leaves = static_cast<long>(leaves.size());
            step.failed_particles = step.leaf_failures > 0;
        }
        stats.merge(step);
        return valid;
    }

   private:
    bool advance(const GuidingCenterPush& pusher, double begin, double end, const GuidingCenter& start,
                 GuidingCenter& finish, int depth, double dt, const IntervalFieldFunction& fields, bool use_aitken,
                 std::vector<Leaf>& leaves, PushStatistics& step) const {
        const double tau = (end - begin) * dt;
        step.depth = std::max(step.depth, depth);
        if (!(tau > 0.0) || !std::isfinite(tau)) {
            if (depth == 0)
                ++step.root_failures;
            return false;
        }
        const auto interval_fields = [&](const Vector3R& a, const Vector3R& b) { return fields(a, b, begin, end); };
        const PushResult result = pusher.push(start, finish, tau, interval_fields, use_aitken);
        step.iteration_sum += result.iterations;
        step.max_iterations = std::max(step.max_iterations, result.iterations);
        if (depth == 0 && !result.converged)
            ++step.root_failures;
        const auto& p = result.state;
        if (!std::isfinite(p.v_parallel) || !std::isfinite(p.position.x()) || !std::isfinite(p.position.y()) ||
            !std::isfinite(p.position.z()) || !std::isfinite(result.position_residual) ||
            !std::isfinite(result.velocity_residual))
            return false;
        finish = p;

        const int limit = enabled_ ? max_depth_ : 0;
        if (!result.converged && depth < limit) {
            const double half = begin + 0.5 * (end - begin);
            if (!(begin < half && half < end))
                return false;
            ++step.retries;
            GuidingCenter mid = start;
            if (!advance(pusher, begin, half, start, mid, depth + 1, dt, fields, use_aitken, leaves, step))
                return false;
            finish = mid;
            return advance(pusher, half, end, mid, finish, depth + 1, dt, fields, use_aitken, leaves, step);
        }
        if (!result.converged) {
            ++step.leaf_failures;
            step.max_leaf_residue_r = std::max(step.max_leaf_residue_r, result.position_residual);
            step.max_leaf_residue_v = std::max(step.max_leaf_residue_v, result.velocity_residual);
        }
        leaves.push_back({begin, end, start, finish});
        return true;
    }

    bool enabled_ = false;
    int max_depth_ = 4;
};

}   // namespace driftKinetic
