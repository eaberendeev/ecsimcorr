#include "boris_pusher.h"
#include "coupling/esirkepov.h"
#include "coupling/gather.h"
#include "simulation_ecsim.h"

namespace algorithmsECSIM {

void predict_velocity_impl_linear(ParticlesArray& particles, const Field3d& fieldEp, const Field3d& fieldB,
                                  const double dt) {
    RECORD_TIMER;
    const double qm = particles.charge / particles.mass();
    const auto& domain = particles.get_domain();
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < particles.size(); ++k) {
        for (auto& particle : particles.particlesData(k)) {
            const Vector3R normalized_coord = domain.to_cell_coordinates(particle.coord);
            const Vector3R E_p = interpolateE_linear(fieldEp, normalized_coord);
            const Vector3R B_p = interpolateB_linear(fieldB, normalized_coord);
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);
        }
    }
}

template <ShapeType ST>
void predict_velocity(ParticlesArray& particles, const Field3d& fieldEp, const Field3d& fieldB, const double dt) {
    RECORD_TIMER;
    if (particles.is_neutral())
        return;

    static_assert(ST == ShapeType::Linear, "Predict velocity supports Linear shape only");
    predict_velocity_impl_linear(particles, fieldEp, fieldB, dt);
}

void predict_current_impl_linear(const ParticlesArray& particles, const Field3d& fieldB, Field3d& fieldJ,
                                 const double dt) {
    RECORD_TIMER;

    const double qp = particles.charge;
    const double mpw = particles.mpw();
    const double q_m = qp / particles.mass();
    const auto& domain = particles.get_domain();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < particles.size(); ++pk) {
        if (particles.particlesData(pk).empty()) {
            continue;
        }
        ParticleShape<Shape, 2> shape, shape05;
        CurrentBuffer<3> cur_buff;
        cur_buff.zero();
        for (auto& particle : particles.particlesData(pk)) {
            const Vector3R cell_coord = domain.to_cell_coordinates(particle.coord);
            shape.fill_from_normalized(cell_coord);
            shape05.fill_from_normalized(cell_coord, Vector3R(0.5, 0.5, 0.5));
            const Vector3R velocity = particle.velocity;

            const Vector3R B_p = interpolateB(fieldB, shape, shape05);

            const Vector3R b = 0.5 * dt * q_m * B_p;

            const double betaI = qp * mpw / (1.0 + b.squaredNorm());

            const Vector3R I_p = betaI * (velocity + velocity.cross(b) + b * velocity.dot(b));

            decompose_current(shape, shape05, I_p, cur_buff);
        }
        auto [start_x, start_y, start_z] = shape.start_.split();
        flush_current_buffer(fieldJ, cur_buff, start_x - 1, start_y - 1, start_z - 1);
    }
}

template <ShapeType ST>
void predict_current(const ParticlesArray& particles, const Field3d& fieldB, Field3d& fieldJ, const double dt) {
    RECORD_TIMER;
    if (particles.is_neutral())
        return;

    static_assert(ST == ShapeType::Linear, "Predict current supports Linear shape only");
    predict_current_impl_linear(particles, fieldB, fieldJ, dt);
}

template void predict_velocity<ShapeType::Linear>(ParticlesArray& particles, const Field3d& fieldEp,
                                                  const Field3d& fieldB, const double dt);
template void predict_current<ShapeType::Linear>(const ParticlesArray& particles, const Field3d& fieldB,
                                                 Field3d& fieldJ, const double dt);

Vector3R calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ, const Grid& grid,
                           const BoundaryConditionHandler& bc) {
    RECORD_TIMER;
    Vector3R potE = Vector3R(0, 0, 0);
    IndexRange range = bc.active_range(grid);

    for (auto i = range.start.x(); i < range.end.x(); ++i) {
        for (auto j = range.start.y(); j < range.end.y(); ++j) {
            for (auto k = range.start.z(); k < range.end.z(); ++k) {
                Vector3R E = Vector3R(fieldE(i, j, k, 0), fieldE(i, j, k, 1), fieldE(i, j, k, 2));
                Vector3R J = Vector3R(fieldJ(i, j, k, 0), fieldJ(i, j, k, 1), fieldJ(i, j, k, 2));
                potE += Vector3R(J.x() * E.x(), J.y() * E.y(), J.z() * E.z());
            }
        }
    }
    return potE;
}

}   // namespace algorithmsECSIM
