#pragma once

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include "particles/boris_pusher.h"
#include "particles/particle_path.h"
#include "whitney_shape.h"

namespace ecsimWhitney {
enum class OrbitPredictor { Ballistic, OldFields, HalfDrift, CenteredFields };

struct FrozenPath {
    Vector3R start, middle, end;   // unwrapped world coordinates
    bool split = false;            // two time intervals, otherwise one chord
    double first_fraction = 0.5;   // split only: fraction of dt on [start, middle]
    double duration = 1.0;         // total fraction of dt covered by the path
    ParticleFate fate = ParticleFate::Unhandled;
    Face wall = Face::UNAVAILABLE;
};

inline FrozenPath predict_path(const Vector3R& start, const Vector3R& velocity, double qm, double dt,
                               const Field3d& old_electric, const Field3d& old_magnetic, const Domain& domain,
                               OrbitPredictor predictor, KernelOrder order = KernelOrder::Linear) {
    const Vector3R middle = start + (0.5 * dt) * velocity;
    if (predictor == OrbitPredictor::Ballistic)
        return {start, middle, start + dt * velocity, false};
    if (predictor == OrbitPredictor::CenteredFields) {
        // v^n is at the CENTER of [x^{n-1/2},x^{n+1/2}]. Predict the
        // velocities at the two quarter-time quadrature points symmetrically.
        // Under uniform E their accelerations cancel in the net displacement.
        const auto e = gather_electric(old_electric, middle, domain, order);
        const auto b = gather_magnetic(old_magnetic, middle, domain, order);
        Particle before(middle, velocity), after(middle, velocity);
        borisPusher::update_vEB(before, qm, e, b, -0.25 * dt);
        borisPusher::update_vEB(after, qm, e, b, 0.25 * dt);
        const auto center = start + (0.5 * dt) * before.velocity;
        return {start, center, center + (0.5 * dt) * after.velocity, true};
    }
    // Only the trial velocity is pushed; the actual ECSIM velocity remains v^n
    // until the field solve. Old fields are held constant in time, not space.
    const Vector3R sample = predictor == OrbitPredictor::HalfDrift ? middle : start;
    Particle trial(sample, velocity);
    borisPusher::update_vEB(trial, qm, gather_electric(old_electric, sample, domain, order),
                            gather_magnetic(old_magnetic, sample, domain, order), dt);
    const auto end = middle + (0.5 * dt) * trial.velocity;
    return {start, middle, end, predictor == OrbitPredictor::HalfDrift};
}

struct OrbitWeights {
    std::vector<EdgeWeight> electric;
    std::vector<EdgeWeight> magnetic;
};

inline void merge_weights(std::vector<EdgeWeight>& weights) {
    std::sort(weights.begin(), weights.end(), [](const auto& a, const auto& b) { return a.index < b.index; });
    size_t count = 0;
    for (const auto& w : weights) {
        if (w.weight == 0)
            continue;
        if (count && weights[count - 1].index == w.index)
            weights[count - 1].weight += w.weight;
        else
            weights[count++] = w;
    }
    weights.resize(count);
}

// Time average of Whitney forms on a straight segment. Split at every grid
// plane (also unwrapped periodic seams). E weights are quadratic within each
// cell, B weights linear; two-point Gauss integration is exact on each piece.
inline void integrate_segment(OrbitWeights& result, const Vector3R& start, const Vector3R& end, double time_fraction,
                              const Domain& domain, bool store_magnetic = true,
                              KernelOrder order = KernelOrder::Linear) {
    const auto u0 = domain.to_cell_coordinates(start), u1 = domain.to_cell_coordinates(end);
    std::vector<double> cuts{0, 1};
    for (int a = 0; a < 3; ++a) {
        if (!std::isfinite(u0[a]) || !std::isfinite(u1[a]))
            throw std::runtime_error("non-finite frozen Whitney trajectory");
        const double lo = std::min(u0[a], u1[a]), hi = std::max(u0[a], u1[a]);
        // Protect integer conversion and prevent unbounded work on a bad orbit.
        if (std::abs(lo) > 1e8 || std::abs(hi) > 1e8 || hi - lo > 1024)
            throw std::runtime_error("frozen Whitney trajectory exceeds integration range");
        for (int plane = static_cast<int>(std::floor(lo)) + 1; plane < hi; ++plane)
            cuts.push_back((plane - u0[a]) / (u1[a] - u0[a]));
    }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    // Linear weights are quadratic in time (2-point Gauss exact); quadratic
    // weights are up to degree four, so use 3-point Gauss.
    static constexpr double gauss2[] = {-1 / std::sqrt(3.0), 1 / std::sqrt(3.0)};
    static constexpr double weight2[] = {1.0, 1.0};
    static constexpr double gauss3[] = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
    static constexpr double weight3[] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    const double* points = order == KernelOrder::Linear ? gauss2 : gauss3;
    const double* quadrature = order == KernelOrder::Linear ? weight2 : weight3;
    const int count = order == KernelOrder::Linear ? 2 : 3;
    for (size_t s = 1; s < cuts.size(); ++s) {
        const double half = 0.5 * (cuts[s] - cuts[s - 1]);
        for (int q = 0; q < count; ++q) {
            const double t = 0.5 * (cuts[s] + cuts[s - 1]) + half * points[q];
            const auto position = start + t * (end - start);
            const double factor = time_fraction * half * quadrature[q];
            const auto electric = electric_stencil(position, domain, order);
            for (int i = 0; i < electric.size; ++i) {
                auto w = electric.weight[i];
                w.weight *= factor;
                result.electric.push_back(w);
            }
            if (store_magnetic) {
                const auto magnetic = magnetic_stencil(position, domain, order);
                for (int i = 0; i < magnetic.size; ++i) {
                    auto w = magnetic.weight[i];
                    w.weight *= factor;
                    result.magnetic.push_back(w);
                }
            }
        }
    }
}

// Clip a one-chord ballistic orbit at the first boundary face selected by the
// PathPolicy. An absorbed path ends on the wall and covers only t_hit of dt.
// A reflected path keeps the incident leg and adds the mirrored outgoing leg;
// the latter is supported for completeness but frozen wall simulations reject
// reflecting conditions until their surface charge is modeled.
inline FrozenPath clip_ballistic_path(const FrozenPath& free, const Particle& particle, const PathPolicy& policy,
                                      const Domain& domain, ParticleFateResult& fate) {
    fate = {};
    if (!policy.active())
        return free;
    Vector3R hit;
    const BoundaryCrossing crossing = shape_particle_path(free.start, free.end, domain, policy, hit);
    if (!crossing.clipped)
        return free;

    Particle at_hit = particle;
    at_hit.coord = crossing.point;
    fate = policy.classify ? policy.classify(at_hit, crossing.face)
                           : ParticleFateResult{ParticleFate::Removed, crossing.face};
    if (fate.fate == ParticleFate::Unhandled)
        fate = {ParticleFate::Removed, crossing.face};

    FrozenPath result = free;
    result.middle = crossing.point;
    result.end = crossing.point;
    result.split = false;
    result.first_fraction = 0.5;
    result.duration = crossing.t_hit;
    result.fate = fate.fate;
    result.wall = crossing.face;
    if (fate.fate == ParticleFate::Reflected) {
        result.end = domain.geom.reflect_from_face(crossing.face, free.end);
        result.split = true;
        result.first_fraction = crossing.t_hit;
        result.duration = 1.0;
    }
    return result;
}

inline OrbitWeights orbit_weights(const FrozenPath& path, const Domain& domain,
                                  KernelOrder order = KernelOrder::Linear) {
    OrbitWeights result;
    const auto add = [&](const Vector3R& start, const Vector3R& end, double fraction) {
        if (fraction > 0)
            integrate_segment(result, start, end, fraction, domain, true, order);
    };
    if (path.split) {
        if (path.first_fraction < 0 || path.first_fraction > 1)
            throw std::invalid_argument("frozen path split fraction must be in [0, 1]");
        add(path.start, path.middle, path.first_fraction);
        add(path.middle, path.end, 1 - path.first_fraction);
    } else {
        if (path.duration < 0 || path.duration > 1)
            throw std::invalid_argument("frozen path duration must be in [0, 1]");
        add(path.start, path.end, path.duration);
    }
    merge_weights(result.electric);
    merge_weights(result.magnetic);
    return result;
}

// One straight part of a full-step particle trajectory. duration is a
// fraction of dt; velocity is the constant physical velocity on this leg.
struct CurrentLeg {
    Vector3R start, end;
    Vector3R velocity;
    double duration = 1.0;
};

// Deposit q/dt * integral(v.W) over the union of straight legs. Each leg is
// integrated separately because the velocity changes at the Boris turn. The
// result is a charge-conserving edge current for CIC density when the legs
// cover the complete particle path from rho^n to rho^{n+1}.
inline void deposit_orbit_current(Field3d& current, double qw, std::initializer_list<CurrentLeg> legs,
                                  const Domain& domain, KernelOrder order = KernelOrder::Linear) {
    // Reuse the per-thread weight buffer: a particle path can cross many grid
    // planes, and allocating one vector per leg dominates this kernel.
    static thread_local OrbitWeights weights;
    weights.magnetic.clear();
    for (const auto& leg : legs) {
        if (leg.duration <= 0)
            continue;
        for (int a = 0; a < 3; ++a)
            if (!std::isfinite(leg.start[a]) || !std::isfinite(leg.end[a]) || !std::isfinite(leg.velocity[a]))
                throw std::runtime_error("non-finite Whitney trajectory current");
        weights.electric.clear();
        integrate_segment(weights, leg.start, leg.end, leg.duration, domain, false, order);
        for (const auto& w : weights.electric) current[w.index] += qw * w.weight * leg.velocity[w.component];
    }
}

// Analytic form of the same orbit integral. Between grid planes the two
// transverse S1 factors are linear in time. For their midpoint values A, B
// and endpoint differences dA, dB, <A(t) B(t)> = A B + dA dB / 12.
// Magnetic S1 weights are linear, so their midpoint values are exact.
inline bool integrate_segment_cell(OrbitWeights& result, const Vector3R& start, const Vector3R& end,
                                   double time_fraction, const Domain& domain, std::vector<double>& cuts) {
    if (time_fraction <= 0) {
        cuts.clear();
        return true;
    }
    const auto u0 = domain.to_cell_coordinates(start), u1 = domain.to_cell_coordinates(end);
    cuts.clear();
    cuts.push_back(0);
    cuts.push_back(1);
    for (int a = 0; a < 3; ++a) {
        if (!std::isfinite(u0[a]) || !std::isfinite(u1[a]))
            throw std::runtime_error("non-finite frozen Whitney trajectory");
        const double lo = std::min(u0[a], u1[a]), hi = std::max(u0[a], u1[a]);
        if (std::abs(lo) > 1e8 || std::abs(hi) > 1e8 || hi - lo > 1024)
            throw std::runtime_error("frozen Whitney trajectory exceeds integration range");
        for (int plane = static_cast<int>(std::floor(lo)) + 1; plane < hi; ++plane)
            cuts.push_back((plane - u0[a]) / (u1[a] - u0[a]));
    }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    const Vector3R displacement = u1 - u0;
    for (size_t s = 1; s < cuts.size(); ++s) {
        const double width = cuts[s] - cuts[s - 1];
        const Vector3R midpoint = u0 + (0.5 * (cuts[s] + cuts[s - 1])) * displacement;
        const Vector3R delta = width * displacement;
        Vector3R fraction;
        Vector3I lower, upper;
        for (int a = 0; a < 3; ++a) {
            int cell = static_cast<int>(std::floor(midpoint[a]));
            fraction[a] = midpoint[a] - cell;
            if (!domain.is_periodic(a)) {
                const int n = domain.num_cells(a);
                if (cell < 0) {
                    cell = 0;
                    fraction[a] = 0;
                } else if (cell >= n) {
                    cell = n - 1;
                    fraction[a] = 1;
                }
            }
            lower[a] = shape_index(cell + domain.grid.ghost_cells(), a, domain);
            upper[a] = shape_index(cell + 1 + domain.grid.ghost_cells(), a, domain);
        }
        const double factor = time_fraction * width;
        for (int c = 0; c < 3; ++c) {
            const int a = (c + 1) % 3, b = (c + 2) % 3;
            for (int ia = 0; ia < 2; ++ia)
                for (int ib = 0; ib < 2; ++ib) {
                    auto node = lower;
                    if (ia)
                        node[a] = upper[a];
                    if (ib)
                        node[b] = upper[b];
                    const double wa = ia ? fraction[a] : 1 - fraction[a];
                    const double wb = ib ? fraction[b] : 1 - fraction[b];
                    const double weight = factor * (wa * wb + ((ia == ib ? 1 : -1) / 12.0) * delta[a] * delta[b]);
                    if (weight != 0)
                        result.electric.push_back({domain.vind(node[0], node[1], node[2], c), c, weight});
                }
            for (int bit = 0; bit < 2; ++bit) {
                auto node = lower;
                if (bit)
                    node[c] = upper[c];
                const double weight = factor * (bit ? fraction[c] : 1 - fraction[c]);
                if (weight != 0)
                    result.magnetic.push_back({domain.vind(node[0], node[1], node[2], c), c, weight});
            }
        }
    }
    return cuts.size() == 2;
}

// Caller-owned scratch is reused between particles. No merge is necessary
// when a straight path stays in one cell: all its edge/face indices differ.
inline void orbit_weights_cell(const FrozenPath& path, const Domain& domain, OrbitWeights& result,
                               std::vector<double>& cuts, KernelOrder order = KernelOrder::Linear) {
    if (order == KernelOrder::Quadratic) {
        // The analytic cell formula below is specialized to two linear factors.
        result = orbit_weights(path, domain, order);
        return;
    }
    result.electric.clear();
    result.magnetic.clear();
    bool one_cell;
    if (path.split) {
        if (path.first_fraction < 0 || path.first_fraction > 1)
            throw std::invalid_argument("frozen path split fraction must be in [0, 1]");
        if (path.first_fraction > 0)
            integrate_segment_cell(result, path.start, path.middle, path.first_fraction, domain, cuts);
        if (path.first_fraction < 1)
            integrate_segment_cell(result, path.middle, path.end, 1 - path.first_fraction, domain, cuts);
        one_cell = false;
    } else {
        if (path.duration < 0 || path.duration > 1)
            throw std::invalid_argument("frozen path duration must be in [0, 1]");
        one_cell = integrate_segment_cell(result, path.start, path.end, path.duration, domain, cuts);
    }
    if (!one_cell) {
        merge_weights(result.electric);
        merge_weights(result.magnetic);
    }
}

inline Vector3R gather_weights(const Field3d& field, const std::vector<EdgeWeight>& weights) {
    Vector3R result(0, 0, 0);
    for (const auto& w : weights) result[w.component] += w.weight * field[w.index];
    return result;
}
}   // namespace ecsimWhitney
