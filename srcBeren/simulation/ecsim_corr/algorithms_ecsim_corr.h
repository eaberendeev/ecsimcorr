#pragma once

#include "ParticlesArray.h"
#include "decompose_esirkepov_current.h"
#include "interpolation.h"

template <ShapeFunction ShapeFn, int ShapeSize>
inline void update_v_move_and_calc_current_impl(ParticlesArray& sp, const Field3d& fieldE, const Field3d& fieldB,
                                                double dt, double half_dt, Field3d& fieldJ, double& pred_w) {
    constexpr auto SMAX = 2 * ShapeSize;

    const double qm = sp.charge / sp.mass();
    const double charge = sp.charge;
    const double mpw = sp.mpw();
    const auto& domain = sp.get_domain();
    auto& particlesData = sp.particlesData;

    const double qx = charge * domain.cell_size().x() / (6 * half_dt) * mpw;
    const double qy = charge * domain.cell_size().y() / (6 * half_dt) * mpw;
    const double qz = charge * domain.cell_size().z() / (6 * half_dt) * mpw;

#pragma omp parallel for schedule(dynamic, 64) reduction(+ : pred_w)
    for (auto pk = 0; pk < particlesData.capacity(); ++pk) {
        auto& particles = particlesData(pk);
        if (particles.empty()) {
            continue;
        }

        ParticleShape<ShapeFn, SMAX> start_shape;
        ParticleShape<ShapeFn, SMAX> end_shape;
        CurrentBuffer<SMAX> curBuf, cellBuf;
        cellBuf.zero();
        curBuf.zero();
        start_shape.fill_zero();

        for (auto& particle : particles) {
            Vector3R start = particle.coord;
            Vector3R old_v = particle.velocity;
            start_shape.fill_from_normalized(domain.to_cell_coordinates(start), GHOST_CELLS);

            const auto norm_coord = domain.to_cell_coordinates(start);
            const Vector3R E_p = interpolateE_linear(fieldE, norm_coord);
            const Vector3R B_p = interpolateB_linear(fieldB, norm_coord);
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);

            particle.move(half_dt);

            Vector3R end = particle.coord;
            end_shape.fill_from_normalized(domain.to_cell_coordinates(end), start_shape.base_, GHOST_CELLS);
            decompose_esirkepov_current(start_shape, end_shape, qx, qy, qz, curBuf);
            cellBuf += curBuf;

            pred_w += mpw * charge * 0.5 * (old_v + particle.velocity).dot(E_p);
        }
        auto [start_x, start_y, start_z] = start_shape.start_.split();
        flush_current_buffer(fieldJ, cellBuf, start_x, start_y, start_z);
    }
}

inline void fused_push_and_deposit(ParticlesArray& sp, const Field3d& fieldE, const Field3d& fieldB, double dt,
                                   double half_dt, Field3d& fieldJ, double& pred_w, ShapeType type) {
    if (sp.is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            std::cout << "fused_push_and_deposit for NGP is not supported\n" << std::endl;
            exit(-1);
        case ShapeType::Linear:
            update_v_move_and_calc_current_impl<Shape, 2>(sp, fieldE, fieldB, dt, half_dt, fieldJ, pred_w);
            break;
        case ShapeType::Quadratic:
            update_v_move_and_calc_current_impl<Shape2, 2>(sp, fieldE, fieldB, dt, half_dt, fieldJ, pred_w);
            break;
    }
}

inline void move_and_calc_current(ParticlesArray& sp, double dt, Field3d& fieldJ, ShapeType type = SHAPE) {
    if (sp.is_neutral())
        return;
    move_and_calc_current(sp.particlesData, sp.get_domain(), sp.charge, sp.mpw(), dt, fieldJ, type);
}

inline void move_and_calc_current_reference(ParticlesArray& sp, double dt, Field3d& fieldJ, ShapeType type = SHAPE) {
    if (sp.is_neutral())
        return;
    move_and_calc_current_reference(sp.particlesData, sp.get_domain(), sp.charge, sp.mpw(), dt, fieldJ, type);
}
