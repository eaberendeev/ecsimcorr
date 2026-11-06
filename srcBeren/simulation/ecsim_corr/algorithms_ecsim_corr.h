#pragma once

#include "coupling/current.h"
#include "coupling/gather.h"
#include "particles/boris_pusher.h"
#include "particles/particle_array.h"

// Accumulate predictor_power = sum(q*mpw*v_mid . E_predictor), without dt.
// fieldJ receives the second half-path current normalized by half_dt; the
// caller averages it with the first half-path current after this kernel.
template <ShapeFunction ShapeFn, int ShapeSize>
inline void update_v_move_and_calc_current_impl(ParticlesArray& sp, const Field3d& fieldE, const Field3d& fieldB,
                                                double dt, double half_dt, Field3d& fieldJ, double& predictor_power,
                                                const PathPolicy& policy) {
    RECORD_TIMER;
    constexpr int SMAX = 2 * ShapeSize;
    static_assert(SMAX == 4, "trajectory deposition uses a common four-node envelope");

    const double qm = sp.charge() / sp.mass();
    const double charge = sp.charge();
    const double mpw = sp.mpw();
    const Domain& domain = sp.get_domain();
    Array3D<std::vector<Particle>>& particlesData = sp.cells();

    const double qx = charge * domain.cell_size().x() / (6 * half_dt) * mpw;
    const double qy = charge * domain.cell_size().y() / (6 * half_dt) * mpw;
    const double qz = charge * domain.cell_size().z() / (6 * half_dt) * mpw;

#pragma omp parallel for schedule(dynamic, 64) reduction(+ : predictor_power)
    for (int pk = 0; pk < particlesData.capacity(); ++pk) {
        std::vector<Particle>& particles = particlesData(pk);
        if (particles.empty()) {
            continue;
        }

        Stencil3D<SMAX> start_shape;
        CurrentBuffer<SMAX> curBuf, cellBuf;
        cellBuf.zero();
        curBuf.zero();

        for (Particle& particle : particles) {
            Vector3R start = particle.coord;
            Vector3R old_v = particle.velocity;
            const Vector3R coord = domain.to_cell_coordinates(start);
            const Vector3I origin(static_cast<int>(std::floor(coord.x())) - 1,
                                  static_cast<int>(std::floor(coord.y())) - 1,
                                  static_cast<int>(std::floor(coord.z())) - 1);
            start_shape = make_stencil_at<ShapeFn, SMAX>(coord, origin);

            const auto no = make_point_stencil<ShapeType::Linear>(coord);
            const auto sh = make_point_stencil<ShapeType::Linear>(coord, Vector3R(0.5));
            const Vector3R E_p = interpolateE(fieldE, no, sh);
            const Vector3R B_p = interpolateB(fieldB, no, sh);
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);

            particle.move(half_dt);

            visit_current_path<ShapeFn>(start_shape, start, particle, domain, policy,
                                        [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool cell_local) {
                                            curBuf.zero();
                                            decompose_esirkepov_current(a, b, qx, qy, qz, curBuf);
                                            if (cell_local)
                                                cellBuf += curBuf;
                                            else
                                                flush_current_buffer(fieldJ, curBuf, a.origin, domain.grid);
                                        });

            predictor_power += mpw * charge * 0.5 * (old_v + particle.velocity).dot(E_p);
        }
        flush_current_buffer(fieldJ, cellBuf, start_shape.origin, domain.grid);
    }
}

template <ShapeType ST = SHAPE>
inline void fused_push_and_deposit(ParticlesArray& sp, const Field3d& fieldE, const Field3d& fieldB, double dt,
                                   double half_dt, Field3d& fieldJ, double& predictor_power,
                                   const PathPolicy& policy = PathPolicy::none()) {
    if (sp.is_neutral())
        return;

    update_v_move_and_calc_current_impl<ShapeKernel<ST>::weight, ShapeKernel<ST>::deposit_half_width>(
        sp, fieldE, fieldB, dt, half_dt, fieldJ, predictor_power, policy);
}
