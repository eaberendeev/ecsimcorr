#pragma once

#include <algorithm>
#include <vector>

#include "coupling/esirkepov.h"
#include "coupling/gather.h"
#include "particles/boris_pusher.h"
#include "particles/particle_array.h"
#include "path_gather.h"
#include "runtime/profiler.h"

// Second half of the ECSIM-Path step: v^n -> v^{n+1} with the mean predicted
// field on the straight step path, x^{n+1/2} -> x^{n+1}. The mean field is
// found by `gather_passes` fixed-point passes from the point field. fieldJ
// receives the straight path current minus the first half-path (already
// deposited by the drift), both normalized by half_dt like the ECSIM-Corr
// second half-path; the caller averages it with the first half. A particle
// whose x^n lies outside the box (a wall changed its first half) keeps the
// point field and the second half-path.
// predictor_power = sum(q*mpw*v_mid . E_used), without dt.
template <ShapeFunction ShapeFn, int ShapeSize>
inline void push_on_path_and_deposit(ParticlesArray& sp, const Field3d& fieldE, const Field3d& fieldB, double dt,
                                     Field3d& fieldJ, double& predictor_power, const PathPolicy& policy,
                                     int gather_passes) {
    RECORD_TIMER;
    constexpr int SMAX = 2 * ShapeSize;
    static_assert(SMAX == 4, "trajectory deposition uses a common four-node envelope");
    if (sp.is_neutral())
        return;

    const double half_dt = 0.5 * dt;
    const double qm = sp.charge() / sp.mass();
    const double charge = sp.charge();
    const double mpw = sp.mpw();
    const Domain& domain = sp.get_domain();
    require_gather_grid(fieldE, domain.grid);
    require_gather_grid(fieldB, domain.grid);
    Array3D<std::vector<Particle>>& particlesData = sp.cells();

    const double qx = charge * domain.cell_size().x() / (6 * half_dt) * mpw;
    const double qy = charge * domain.cell_size().y() / (6 * half_dt) * mpw;
    const double qz = charge * domain.cell_size().z() / (6 * half_dt) * mpw;
    // Trial paths and the removed first half-path must not record wall events.
    PathPolicy trial_policy = policy;
    trial_policy.events = nullptr;
    bool periodic[3] = {policy.periodic[0], policy.periodic[1], policy.periodic[2]};
    const Vector3R cell = domain.cell_size();
    const double tiny = 1e-8 * std::min(cell.x(), std::min(cell.y(), cell.z()));

#pragma omp parallel for schedule(dynamic, 64) reduction(+ : predictor_power)
    for (int pk = 0; pk < particlesData.capacity(); ++pk) {
        std::vector<Particle>& particles = particlesData(pk);
        if (particles.empty())
            continue;

        const Vector3I point_origin = checked_floor_index(domain.to_cell_coordinates(particles.front().coord));
        const Vector3I cell_origin = deposition_origin(point_origin);
        CurrentBuffer<SMAX> cellBuf;
        cellBuf.zero();
        std::vector<pathGather::WeightEntry> weights;

        for (Particle& particle : particles) {
            const Vector3R start = particle.coord;
            const Vector3R old_v = particle.velocity;
            const Vector3R coord = domain.to_cell_coordinates(start);
            const Stencil3D<2> no = make_linear_stencil_at(coord, point_origin);
            const Stencil3D<2> sh =
                make_linear_stencil_at(coord, linear_half_shifted_origin(coord, point_origin), Vector3R::Constant(0.5));
            Vector3R E_p = interpolateE(fieldE, no, sh);
            const Vector3R B_p = interpolateB(fieldB, no, sh);

            Vector3R path_start = start - half_dt * old_v;
            Vector3R first_end = start;
            pathGather::unwrap_into_box(domain, policy, path_start, first_end);
            const bool straight = domain.geom.contains(path_start, periodic);
            if (straight)
                for (int pass = 0; pass < gather_passes; ++pass) {
                    Particle trial = particle;
                    borisPusher::update_vEB(trial, qm, E_p, B_p, dt);
                    pathGather::path_weights<ShapeFn>(path_start, half_dt * (old_v + trial.velocity), domain,
                                                      trial_policy, dt, tiny, weights);
                    E_p = pathGather::gather(fieldE, weights);
                }
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);
            particle.move(half_dt);

            if (straight) {
                // Add the straight step path, take the first half-path away
                // (both conserve charge).
                const Vector3R end = particle.coord;
                particle.coord = first_end + (end - start);
                pathGather::deposit_path<ShapeFn>(fieldJ, path_start, particle, domain, policy, qx, qy, qz);
                particle.coord = end;
                Particle first = particle;
                first.coord = first_end;
                first.velocity = old_v;
                pathGather::deposit_path<ShapeFn>(fieldJ, path_start, first, domain, trial_policy, -qx, -qy, -qz);
            } else {
                const Stencil3D<SMAX> start_shape = make_stencil_at<ShapeFn, SMAX>(coord, cell_origin);
                visit_current_path<ShapeFn>(
                    start_shape, start, particle, domain, policy,
                    [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool cell_local) {
                        if (cell_local) {
                            decompose_esirkepov_current(a, b, qx, qy, qz, cellBuf);
                        } else {
                            CurrentBuffer<SMAX> scratch;
                            scratch.zero();
                            decompose_esirkepov_current(a, b, qx, qy, qz, scratch);
                            flush_current_buffer(fieldJ, scratch, to_storage_index(a.origin));
                        }
                    });
            }
            predictor_power += mpw * charge * 0.5 * (old_v + particle.velocity).dot(E_p);
        }
        flush_current_buffer(fieldJ, cellBuf, to_storage_index(cell_origin));
    }
}
