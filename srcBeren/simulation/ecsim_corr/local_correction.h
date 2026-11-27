#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

#include "coupling/esirkepov.h"
#include "coupling/gather.h"
#include "coupling/stencil.h"
#include "fields/field_boundaries.h"
#include "particles/boris_pusher.h"
#include "particles/particle_array.h"
#include "particles/particle_path.h"
#include "runtime/profiler.h"

// Per-particle energy correction of ECSIM-Corr. Every particle receives the
// work of the corrected field (E^n + E^{n+1})/2 + E_ext on its own Esirkepov
// path current, minus its predictor work q dt v_mid . E_p, as a kick along the
// force the predictor missed (E_corr - E_p at x^{n+1/2}). The sum of these
// defects is the global correction target, so a species factor is needed only
// for the remainder (paths reflected in the first half, kicks that cannot be
// made, work at wall nodes masked out of the current).
//
// Nothing is stored per particle: x^{n+1/2} = x^{n+1} - dt/2 v^{n+1}, E_p and
// B_p are gathered there again with the predictor stencils, v^n follows from
// the time-reversible velocity update with -E_p, -B_p, and x^n = x^{n+1/2} -
// dt/2 v^n. Both half-paths follow the same visit_current_path segments as the deposit.
namespace ecsimCorrLocal {

struct Result {
    double local = 0;     // energy given by per-particle kicks
    int64_t pooled = 0;   // particles left to the species factor
    double kinetic = 0;   // kinetic energy of the species after the kicks
};

// The particles in a bin share the midpoint cell. Pack its corrected field
// once, with each component arranged in three planes normal to that axis.
// Nodes outside storage contribute zero, as in flush_current_buffer.
struct CellWorkField {
    Vector3I origin;
    alignas(64) double planes[3][3][4][4];

    CellWorkField(const Field3d& field, const Vector3I& stencil_origin) : origin(stencil_origin) {
        const Vector3I base = to_storage_index(origin);
        const auto sample = [&](int x, int y, int z, Axis axis) {
            const int i = base.x() + x, j = base.y() + y, k = base.z() + z;
            return field.checkIndex(i, j, k, static_cast<int>(axis)) ? field(i, j, k, axis) : 0.0;
        };
        for (int plane = 0; plane < 3; ++plane)
            for (int u = 0; u < 4; ++u)
                for (int v = 0; v < 4; ++v) {
                    planes[0][plane][u][v] = sample(plane, u, v, Axis::X);
                    planes[1][plane][u][v] = sample(v, plane, u, Axis::Y);
                    planes[2][plane][u][v] = sample(u, v, plane, Axis::Z);
                }
    }

    // The Esirkepov current along each axis is a prefix sum of the change in
    // shape times a transverse weight. Apply those weights directly to E.
    double work(const Stencil3D<4>& start, const Stencil3D<4>& end, const Vector3R& charge_factors) const {
        double total = 0;
        for (int d = 0; d < 3; ++d) {
            const Axis along = static_cast<Axis>(d);
            const Axis across_u = static_cast<Axis>((d + 1) % 3);
            const Axis across_v = static_cast<Axis>((d + 2) % 3);
            const double p0 = start(0, along) - end(0, along);
            const double p1 = p0 + (start(1, along) - end(1, along));
            const double p2 = p1 + (start(2, along) - end(2, along));
            const auto& e = planes[d];
            double sums[4] = {};
            for (int u = 0; u < 4; ++u) {
                const double start_weight = 2 * start(u, across_u) + end(u, across_u);
                const double end_weight = 2 * end(u, across_u) + start(u, across_u);
#pragma omp simd
                for (int v = 0; v < 4; ++v) {
                    const double transverse = start_weight * start(v, across_v) + end_weight * end(v, across_v);
                    sums[v] += transverse * (p0 * e[0][u][v] + p1 * e[1][u][v] + p2 * e[2][u][v]);
                }
            }
            total += charge_factors[d] * (sums[0] + sums[1] + sums[2] + sums[3]);
        }
        return total;
    }
};

// Work of `field` on the Esirkepov current of the straight path start -> end
// (normalized by `half_dt`, as in the deposit), including wall clipping.
template <ShapeFunction ShapeFn>
inline double path_work(const Field3d& field, const Vector3R& start, const Particle& trial, const Domain& domain,
                        const PathPolicy& policy, const Vector3R& charge_factors,
                        const CellWorkField* cell_field = nullptr) {
    constexpr int SMAX = 4;
    const Stencil3D<SMAX> start_shape = make_deposition_stencil<ShapeFn, SMAX>(domain.to_cell_coordinates(start));
    double work = 0;
    visit_current_path<ShapeFn>(start_shape, start, trial, domain, policy,
                                [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool) {
                                    if (cell_field && cell_field->origin == a.origin) {
                                        work += cell_field->work(a, b, charge_factors);
                                        return;
                                    }
                                    work += CellWorkField(field, a.origin).work(a, b, charge_factors);
                                });
    return work;
}

// A quadratic shape that fits a four-node stencil is centered on one of
// its two inner nodes (ties at the outer limits have a zero end weight).
template <ShapeFunction ShapeFn>
inline Stencil3D<4> fitted_shape(const Vector3R& coord, const Vector3I& origin) {
    if constexpr (ShapeFn != Shape2)
        return make_stencil_at<ShapeFn, 4>(coord, origin);
    Stencil3D<4> shape;
    shape.origin = origin;
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        const int lower = origin[axis] + 1;
        const bool upper = coord[axis] > lower + 0.5;
        const double delta = coord[axis] - (upper ? lower + 1 : lower);
        const double left = 0.5 * (0.5 - delta) * (0.5 - delta);
        const double center = 0.75 - delta * delta;
        const double right = 0.5 * (0.5 + delta) * (0.5 + delta);
        shape(0, axis) = upper ? 0.0 : left;
        shape(1, axis) = upper ? left : center;
        shape(2, axis) = upper ? center : right;
        shape(3, axis) = upper ? right : 0.0;
    }
    return shape;
}

// For short straight paths all three shapes can share the midpoint field.
template <ShapeFunction ShapeFn>
inline double two_path_work(const Field3d& field, const CellWorkField& cell, const Domain& domain,
                            const PathPolicy& policy, const Vector3R& charge_factors, const Vector3R& start,
                            const Particle& first, const Vector3R& middle, const Vector3R& mid_coord,
                            const Particle& last) {
    const auto mid_shape = fitted_shape<ShapeFn>(mid_coord, deposition_origin(mid_coord));
    const Vector3R begin = domain.to_cell_coordinates(start);
    const Vector3R end = domain.to_cell_coordinates(last.coord);
    // The first deposit may have been anchored in a neighboring cell. At
    // most one cell of motion guarantees it didn't split the segment there.
    if (first.coord == middle && mid_shape.origin == cell.origin && (mid_coord - begin).abs_max() <= 1.0 &&
        stencil_support_fits<ShapeFn>(mid_shape, begin) && stencil_support_fits<ShapeFn>(mid_shape, end)) {
        Vector3R clipped;
        if (!policy.active() || (!shape_particle_path(start, middle, domain, policy, clipped).clipped &&
                                 !shape_particle_path(middle, last.coord, domain, policy, clipped).clipped)) {
            const auto a = fitted_shape<ShapeFn>(begin, cell.origin);
            const auto b = fitted_shape<ShapeFn>(end, cell.origin);
            return cell.work(a, mid_shape, charge_factors) + cell.work(mid_shape, b, charge_factors);
        }
    }
    const double first_work = path_work<ShapeFn>(field, start, first, domain, policy, charge_factors, &cell);
    return first_work + path_work<ShapeFn>(field, middle, last, domain, policy, charge_factors, &cell);
}

// Shift a periodic path so that its start lies in the box, as in the first
// half-drift (which started from x^n inside the box).
inline void unwrap_into_box(const Domain& domain, const FieldBoundaryConditions& boundaries, Vector3R& start,
                            Vector3R& end) {
    const Vector3R lower = domain.geom.box_min, upper = domain.geom.box_max;
    for (int a = 0; a < 3; ++a) {
        if (!boundaries.is_periodic(a))
            continue;
        const double length = upper[a] - lower[a];
        while (start[a] < lower[a]) {
            start[a] += length;
            end[a] += length;
        }
        while (start[a] >= upper[a]) {
            start[a] -= length;
            end[a] -= length;
        }
    }
}

template <ShapeFunction ShapeFn>
inline Result correct(ParticlesArray& sp, const Field3d& predictor_electric, const Field3d& magnetic,
                      const Field3d& corrected_electric, double dt, const PathPolicy& policy,
                      const FieldBoundaryConditions& boundaries) {
    RECORD_TIMER;
    Result total;
    if (sp.is_neutral())
        return total;
    const Domain& domain = sp.get_domain();
    const double qm = sp.charge() / sp.mass();
    const double qw = sp.charge() * sp.mpw();
    const double mw = sp.mass() * sp.mpw();
    const double half_dt = 0.5 * dt;
    // Deposit normalization of one half-path; the step current is their mean.
    const Vector3R charge_factors = sp.charge() * domain.cell_size() / (6 * half_dt) * sp.mpw();
    const bool periodic[3] = {boundaries.is_periodic(0), boundaries.is_periodic(1), boundaries.is_periodic(2)};
    double local = 0;
    int64_t pooled = 0;
    double kinetic = 0;
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : local, pooled, kinetic)
    for (int cell = 0; cell < sp.size(); ++cell) {
        if (sp(cell).empty())
            continue;
        const Vector3I bin(domain.pos_sind(cell, Axis::X), domain.pos_sind(cell, Axis::Y),
                           domain.pos_sind(cell, Axis::Z));
        const CellWorkField cell_field(corrected_electric, deposition_origin(to_node_index(bin)));
        for (Particle& p : sp(cell)) {
            const Vector3R middle = p.coord - half_dt * p.velocity;
            const Vector3R coord = domain.to_cell_coordinates(middle);
            const Vector3I origin = checked_floor_index(coord);
            const Stencil3D<2> no = make_linear_stencil_at(coord, origin);
            const Stencil3D<2> sh =
                make_linear_stencil_at(coord, linear_half_shifted_origin(coord, origin), Vector3R::Constant(0.5));
            const Vector3R e_p = interpolateE(predictor_electric, no, sh);
            const Vector3R b_p = interpolateB(magnetic, no, sh);
            Particle first = p;
            borisPusher::update_vEB(first, qm, -1.0 * e_p, -1.0 * b_p, dt);
            const Vector3R old_velocity = first.velocity;
            // First half-path x^n -> x^{n+1/2} with v^n. A start outside the
            // domain means a wall changed the first half: leave it pooled.
            Vector3R start = middle - half_dt * old_velocity, mid = middle;
            unwrap_into_box(domain, boundaries, start, mid);
            if (!domain.geom.contains(start, periodic)) {
                ++pooled;
                continue;
            }
            first.coord = mid;
            // Absorbed particles receive their kick before the boundary pass removes them.
            const double corrected_work =
                half_dt * two_path_work<ShapeFn>(corrected_electric, cell_field, domain, policy, charge_factors, start,
                                                 first, middle, coord, p);
            const double predictor_work = dt * qw * 0.5 * (old_velocity + p.velocity).dot(e_p);
            const double defect = corrected_work - predictor_work;
            // Kick along the force the predictor missed, else along v.
            Vector3R direction = interpolateE(corrected_electric, no, sh) - e_p;
            if (direction.norm() <= 1e-14 * (e_p.norm() + std::numeric_limits<double>::min()))
                direction = p.velocity;
            const double norm = direction.norm();
            if (!(norm > 0)) {
                ++pooled;
                continue;
            }
            direction = direction / norm;
            double along = p.velocity.dot(direction);
            double discriminant = along * along + 2 * defect / mw;
            if (discriminant < 0) {
                const double speed = p.velocity.norm();
                if (!(speed > 0)) {
                    ++pooled;
                    continue;
                }
                direction = p.velocity / speed;
                along = speed;
                discriminant = along * along + 2 * defect / mw;
                if (discriminant < 0) {
                    // Stop it; the rest goes to the species factor.
                    local -= 0.5 * mw * speed * speed;
                    p.velocity = Vector3R(0, 0, 0);
                    ++pooled;
                    continue;
                }
            }
            const double root = std::sqrt(discriminant);
            p.velocity += (along >= 0 ? root - along : -along - root) * direction;
            local += defect;
        }
        // The particles of the cell are still in cache.
        for (const Particle& p : sp(cell)) kinetic += get_energy_particle(p.velocity, sp.mass(), sp.mpw());
    }
    total.local = local;
    total.pooled = pooled;
    total.kinetic = kinetic;
    return total;
}

}   // namespace ecsimCorrLocal
