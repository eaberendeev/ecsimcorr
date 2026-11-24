#pragma once

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "coupling/esirkepov.h"
#include "coupling/gather.h"
#include "particles/particle_array.h"
#include "particles/particle_path.h"
#include "runtime/profiler.h"

// Fields and currents on the straight step path of ECSIM-Path (after Chen,
// Chacon, Barnes, JCP 230 (2011)). A particle is pushed with the mean field on
// x^n -> x^{n+1} = x^n + dt v_mid and deposits that path, so its energy changes
// by the work of the field on its own current: q w (d_c / dt) E_c summed over
// the Esirkepov weights w of the path, d = dt v_mid. (Two half-paths would not
// do: under gyration they go out and back, d_c -> 0 while their works do not
// cancel.)
namespace pathGather {

// Shift a path so that its start lies in the box along the periodic axes, as
// the first half-drift started from x^n inside the box.
inline void unwrap_into_box(const Domain& domain, const PathPolicy& policy, Vector3R& start, Vector3R& end) {
    const Vector3R lower = domain.geom.box_min, upper = domain.geom.box_max;
    for (int axis = 0; axis < 3; ++axis) {
        if (!policy.periodic[axis])
            continue;
        const double length = upper[axis] - lower[axis];
        const double periods = std::floor((start[axis] - lower[axis]) / length);
        start[axis] -= periods * length;
        end[axis] -= periods * length;
    }
}

struct WeightEntry {
    Vector3I node;   // storage index
    int component;
    double weight;
};

// Esirkepov weights of the straight path start -> start + d by component: the
// step current of a particle q w moving on it is q w (d_c / dt) sum w at the
// nodes, so sum w E_c(node) is the mean field component c on the path (the
// weights of one component sum to 1). Only the components marked in `take`.
template <ShapeFunction ShapeFn>
inline void append_path_weights(const Vector3R& start, const Vector3R& d, const bool take[3], const Domain& domain,
                                const PathPolicy& policy, double dt, std::vector<WeightEntry>& out) {
    constexpr int SMAX = 4;
    const Vector3R cell = domain.cell_size();
    Particle trial;
    trial.coord = start + d;
    trial.velocity = d / dt;
    const Stencil3D<SMAX> start_shape = make_deposition_stencil<ShapeFn, SMAX>(domain.to_cell_coordinates(start));
    visit_current_path<ShapeFn>(
        start_shape, start, trial, domain, policy, [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool) {
            CurrentBuffer<SMAX> scratch;
            scratch.zero();
            decompose_esirkepov_current(a, b, cell.x() / (6 * dt), cell.y() / (6 * dt), cell.z() / (6 * dt), scratch);
            const Vector3I origin = to_storage_index(a.origin);
            for (int n = 0; n < SMAX; ++n)
                for (int m = 0; m < SMAX; ++m)
                    for (int k = 0; k < SMAX; ++k)
                        for (int c = 0; c < 3; ++c) {
                            const double current = scratch(n, m, k, c);
                            if (take[c] && current != 0)
                                out.push_back({origin + Vector3I(n, m, k), c, current * dt / d[c]});
                        }
        });
}

// A component the path hardly moves along (|d_c| <= tiny) takes the weights of
// the path moved by tiny along it: the weights are Lipschitz in d_c, so this is
// their limit to O(tiny). The policy must not record wall events.
template <ShapeFunction ShapeFn>
inline void path_weights(const Vector3R& start, const Vector3R& d, const Domain& domain, const PathPolicy& policy,
                         double dt, double tiny, std::vector<WeightEntry>& out) {
    out.clear();
    const bool regular[3] = {std::abs(d.x()) > tiny, std::abs(d.y()) > tiny, std::abs(d.z()) > tiny};
    if (regular[0] || regular[1] || regular[2])
        append_path_weights<ShapeFn>(start, d, regular, domain, policy, dt, out);
    for (int c = 0; c < 3; ++c) {
        if (regular[c])
            continue;
        Vector3R moved = d;
        moved[c] = tiny;
        const bool take[3] = {c == 0, c == 1, c == 2};
        append_path_weights<ShapeFn>(start, moved, take, domain, policy, dt, out);
    }
}

inline bool inside_storage(const Vector3I& node, const Vector3I& sizes) {
    return node.x() >= 0 && node.y() >= 0 && node.z() >= 0 && node.x() < sizes.x() && node.y() < sizes.y() &&
           node.z() < sizes.z();
}

// Mean field on the path: sum w E_c(node); nodes outside the storage are zero
// padding, as in flush_current_buffer.
inline Vector3R gather(const Field3d& field, const std::vector<WeightEntry>& weights) {
    Vector3R mean = Vector3R::Zero();
    const Vector3I sizes = field.sizes();
    for (const WeightEntry& e : weights)
        if (inside_storage(e.node, sizes))
            mean[e.component] += e.weight * field(e.node.x(), e.node.y(), e.node.z(), e.component);
    return mean;
}

// Deposit the Esirkepov current of start -> trial.coord with the given
// normalization (negative to take a path away), scattered atomically.
template <ShapeFunction ShapeFn>
inline void deposit_path(Field3d& fieldJ, const Vector3R& start, const Particle& trial, const Domain& domain,
                         const PathPolicy& policy, double qx, double qy, double qz) {
    constexpr int SMAX = 4;
    const Stencil3D<SMAX> start_shape = make_deposition_stencil<ShapeFn, SMAX>(domain.to_cell_coordinates(start));
    visit_current_path<ShapeFn>(start_shape, start, trial, domain, policy,
                                [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool) {
                                    CurrentBuffer<SMAX> scratch;
                                    scratch.zero();
                                    decompose_esirkepov_current(a, b, qx, qy, qz, scratch);
                                    flush_current_buffer(fieldJ, scratch, to_storage_index(a.origin));
                                });
}

// Predictor current of the path that the push will deposit: the weights are
// frozen on the path predicted without the field, d = dt alpha v^n from x^n =
// x^{n+1/2} - dt/2 v^n, J0 = q w w_c (alpha v^n)_c. The ECSIM response L stays
// the point one (a path response gives the same energy error at a much higher
// cost, both freeze the weights).
template <ShapeFunction ShapeFn>
inline void predict_current(const ParticlesArray& sp, const Field3d& fieldB, const PathPolicy& policy, double dt,
                            Field3d& fieldJ) {
    RECORD_TIMER;
    if (sp.is_neutral())
        return;
    const Domain& domain = sp.get_domain();
    const double qm = sp.charge() / sp.mass();
    const double qw = sp.charge() * sp.mpw();
    const double half_dt = 0.5 * dt;
    const Vector3R cell = domain.cell_size();
    const double tiny = 1e-8 * std::min(cell.x(), std::min(cell.y(), cell.z()));
    bool periodic[3] = {policy.periodic[0], policy.periodic[1], policy.periodic[2]};
    PathPolicy trial_policy = policy;
    trial_policy.events = nullptr;
    const Vector3I sizes = fieldJ.sizes();
#pragma omp parallel
    {
        std::vector<WeightEntry> weights;
#pragma omp for schedule(dynamic, sp.omp_granularity())
        for (int pk = 0; pk < sp.size(); ++pk)
            for (const Particle& p : sp(pk)) {
                const Vector3R coord = domain.to_cell_coordinates(p.coord);
                const Vector3I origin = checked_floor_index(coord);
                const Stencil3D<2> no = make_linear_stencil_at(coord, origin);
                const Stencil3D<2> sh = make_linear_stencil_at(coord, linear_half_shifted_origin(coord, origin),
                                                               Vector3R::Constant(0.5));
                const Vector3R b = half_dt * qm * interpolateB(fieldB, no, sh);
                const Vector3R& v = p.velocity;
                const Vector3R alpha_v = (v + v.cross(b) + b * v.dot(b)) / (1.0 + b.squaredNorm());
                const Vector3R d = dt * alpha_v;
                Vector3R start = p.coord - half_dt * v, end = p.coord;
                unwrap_into_box(domain, policy, start, end);
                if (!domain.geom.contains(start, periodic))
                    start = p.coord - 0.5 * d;
                path_weights<ShapeFn>(start, d, domain, trial_policy, dt, tiny, weights);
                for (const WeightEntry& e : weights)
                    if (inside_storage(e.node, sizes)) {
#pragma omp atomic
                        fieldJ(e.node.x(), e.node.y(), e.node.z(), e.component) += qw * e.weight * alpha_v[e.component];
                    }
            }
    }
}

}   // namespace pathGather
