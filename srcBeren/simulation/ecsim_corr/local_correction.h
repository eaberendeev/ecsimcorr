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
};

// The particles in a bin share the midpoint cell. Pack its corrected field
// once, with each component arranged in three planes normal to that axis.
// Nodes outside storage contribute zero, as in flush_current_buffer.
struct CellWorkField {
    Vector3I origin;
    alignas(64) double planes[3][3][4][4];

    CellWorkField(const Field3d& field, const Vector3I& stencil_origin) : origin(stencil_origin) {
        const Vector3I base = to_storage_index(origin);
        const Vector3I sizes = field.sizes();
        for (int x = 0; x < 4; ++x)
            for (int y = 0; y < 4; ++y)
                for (int z = 0; z < 4; ++z) {
                    const int i = base.x() + x, j = base.y() + y, k = base.z() + z;
                    const bool inside = i >= 0 && i < sizes.x() && j >= 0 && j < sizes.y() && k >= 0 && k < sizes.z();
                    if (x < 3)
                        planes[0][x][y][z] = inside ? field(i, j, k, 0) : 0;
                    if (y < 3)
                        planes[1][y][z][x] = inside ? field(i, j, k, 1) : 0;
                    if (z < 3)
                        planes[2][z][x][y] = inside ? field(i, j, k, 2) : 0;
                }
    }

    // In the Esirkepov recurrence the transverse factor is constant down
    // a column. Contract the three longitudinal prefix sums with E first,
    // then apply the transverse factor directly to that scalar work.
    template <Axis D, Axis U, Axis V>
    double component(const Stencil3D<4>& a, const Stencil3D<4>& b, double q) const {
        const double p0 = a(0, D) - b(0, D);
        const double p1 = p0 + (a(1, D) - b(1, D));
        const double p2 = p1 + (a(2, D) - b(2, D));
        const auto& e = planes[static_cast<int>(D)];
        double sum = 0;
        for (int u = 0; u < 4; ++u) {
            const double at = 2 * a(u, U) + b(u, U);
            const double bt = 2 * b(u, U) + a(u, U);
#pragma omp simd reduction(+ : sum)
            for (int v = 0; v < 4; ++v) {
                const double transverse = at * a(v, V) + bt * b(v, V);
                sum += transverse * (p0 * e[0][u][v] + p1 * e[1][u][v] + p2 * e[2][u][v]);
            }
        }
        return q * sum;
    }

    double work(const Stencil3D<4>& a, const Stencil3D<4>& b, double qx, double qy, double qz) const {
        return component<Axis::X, Axis::Y, Axis::Z>(a, b, qx) + component<Axis::Y, Axis::Z, Axis::X>(a, b, qy) +
               component<Axis::Z, Axis::X, Axis::Y>(a, b, qz);
    }
};

// Work of `field` on the Esirkepov current of the straight path start -> end
// (normalized by `half_dt`, as in the deposit), including wall clipping.
template <ShapeFunction ShapeFn>
inline double path_work(const Field3d& field, const Vector3R& start, const Particle& trial, const Domain& domain,
                        const PathPolicy& policy, double qx, double qy, double qz,
                        const CellWorkField* cell_field = nullptr) {
    constexpr int SMAX = 4;
    const Stencil3D<SMAX> start_shape = make_deposition_stencil<ShapeFn, SMAX>(domain.to_cell_coordinates(start));
    double work = 0;
    visit_current_path<ShapeFn>(start_shape, start, trial, domain, policy,
                                [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool) {
                                    if (cell_field && cell_field->origin == a.origin) {
                                        work += cell_field->work(a, b, qx, qy, qz);
                                        return;
                                    }
                                    work += CellWorkField(field, a.origin).work(a, b, qx, qy, qz);
                                });
    return work;
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
    const double qx = sp.charge() * domain.cell_size().x() / (6 * half_dt) * sp.mpw();
    const double qy = sp.charge() * domain.cell_size().y() / (6 * half_dt) * sp.mpw();
    const double qz = sp.charge() * domain.cell_size().z() / (6 * half_dt) * sp.mpw();
    const bool periodic[3] = {boundaries.is_periodic(0), boundaries.is_periodic(1), boundaries.is_periodic(2)};
    double local = 0;
    int64_t pooled = 0;
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : local, pooled)
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
            first.velocity = old_velocity;
            const double work_first =
                path_work<ShapeFn>(corrected_electric, start, first, domain, policy, qx, qy, qz, &cell_field);
            const double work_second =
                path_work<ShapeFn>(corrected_electric, middle, p, domain, policy, qx, qy, qz, &cell_field);
            // Absorbed particles are corrected too: the boundary pass removes
            // them afterwards with their corrected kinetic energy.
            const double corrected_work = dt * 0.5 * (work_first + work_second);
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
    }
    total.local = local;
    total.pooled = pooled;
    return total;
}

}   // namespace ecsimCorrLocal
