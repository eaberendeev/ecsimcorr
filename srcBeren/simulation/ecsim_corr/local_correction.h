#pragma once

#include <cmath>
#include <limits>

#include "coupling/esirkepov.h"
#include "coupling/gather.h"
#include "coupling/stencil.h"
#include "fields/field_boundaries.h"
#include "particles/boris_pusher.h"
#include "particles/particle_array.h"
#include "particles/particle_path.h"

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
// dt/2 v^n. The Esirkepov buffers of both half-paths are rebuilt with the same
// visit_current_path/decompose_esirkepov_current as the deposit.
namespace ecsimCorrLocal {

struct Result {
    double local = 0;       // energy given by per-particle kicks
    long long pooled = 0;   // particles left to the species factor
};

// sum over the buffer of J . E at the storage nodes the deposit would reach.
template <int SMAX>
inline double buffer_work(const Field3d& field, const CurrentBuffer<SMAX>& buf, const Vector3I& storage_origin) {
    double work = 0;
    const Vector3I sizes = field.sizes();
    for (int n = 0; n < SMAX; ++n) {
        const int ix = storage_origin.x() + n;
        if (ix < 0 || ix >= sizes.x())
            continue;
        for (int m = 0; m < SMAX; ++m) {
            const int iy = storage_origin.y() + m;
            if (iy < 0 || iy >= sizes.y())
                continue;
            for (int k = 0; k < SMAX; ++k) {
                const int iz = storage_origin.z() + k;
                if (iz < 0 || iz >= sizes.z())
                    continue;
                for (int d = 0; d < 3; ++d) work += buf(n, m, k, d) * field(ix, iy, iz, d);
            }
        }
    }
    return work;
}

// Work of `field` on the Esirkepov current of the straight path start -> end
// (normalized by `half_dt`, as in the deposit), including wall clipping.
template <ShapeFunction ShapeFn>
inline double path_work(const Field3d& field, const Vector3R& start, const Particle& trial, const Domain& domain,
                        const PathPolicy& policy, double qx, double qy, double qz) {
    constexpr int SMAX = 4;
    const Stencil3D<SMAX> start_shape = make_deposition_stencil<ShapeFn, SMAX>(domain.to_cell_coordinates(start));
    double work = 0;
    visit_current_path<ShapeFn>(start_shape, start, trial, domain, policy,
                                [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool) {
                                    CurrentBuffer<SMAX> scratch;
                                    scratch.zero();
                                    decompose_esirkepov_current(a, b, qx, qy, qz, scratch);
                                    work += buffer_work(field, scratch, to_storage_index(a.origin));
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
    bool periodic[3] = {boundaries.is_periodic(0), boundaries.is_periodic(1), boundaries.is_periodic(2)};
    double local = 0;
    long long pooled = 0;
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : local, pooled)
    for (int cell = 0; cell < sp.size(); ++cell) {
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
            const double work_first = path_work<ShapeFn>(corrected_electric, start, first, domain, policy, qx, qy, qz);
            const double work_second = path_work<ShapeFn>(corrected_electric, middle, p, domain, policy, qx, qy, qz);
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
