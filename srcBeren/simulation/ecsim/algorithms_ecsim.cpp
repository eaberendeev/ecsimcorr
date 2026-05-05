#include "interpolation.h"
#include "simulation_ecsim.h"

namespace algorithmsECSIM {

void predict_velocity_impl_linear(ParticlesArray& particles,
                                  const Field3d& fieldEp, const Field3d& fieldB,
                                  const double dt) {
    const double qm = particles.charge / particles.mass();
    const auto& domain = particles.get_domain();
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < particles.size(); ++k) {
        for (auto& particle : particles.particlesData(k)) {
            const auto normalized_coord =
                domain.to_cell_coordinates(particle.coord);
            const Vector3R E_p = interpolateE_linear(fieldEp, normalized_coord);
            const Vector3R B_p = interpolateB_linear(fieldB, normalized_coord);
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);
        }
    }
}

void predict_velocity_impl_ngp(ParticlesArray& particles,
                               const Field3d& fieldEp, const Field3d& fieldB,
                               const double dt) {
    const double qm = particles.charge / particles.mass();
    const auto& domain = particles.get_domain();

#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < particles.size(); ++k) {
        for (auto& particle : particles.particlesData(k)) {
            const auto normalized_coord =
                domain.to_cell_coordinates(particle.coord);
            const Vector3R E_p = interpolateE_ngp(fieldEp, normalized_coord);
            const Vector3R B_p = interpolateB_ngp(fieldB, normalized_coord);
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);
        }
    }
}

void predict_velocity(ParticlesArray& particles, const Field3d& fieldEp,
                      const Field3d& fieldB, const double dt, ShapeType type) {
    if (particles.is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            predict_velocity_impl_ngp(particles, fieldEp, fieldB, dt);
            break;
        case ShapeType::Linear:
            predict_velocity_impl_linear(particles, fieldEp, fieldB, dt);
            break;
        case ShapeType::Quadratic:
            std::cout
                << "Predict velocity for quadratic shape not implemented\n";
            exit(-1);
    }
}

void predict_current_impl_linear(const ParticlesArray& particles,
                                 const Field3d& fieldB, Field3d& fieldJ,
                                 const double dt) {
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
            const Vector3R cell_coord =
                domain.to_cell_coordinates(particle.coord);
            shape.fill_from_normalized(cell_coord);
            shape05.fill_from_normalized(cell_coord, Vector3R(0.5, 0.5, 0.5));
            const Vector3R velocity = particle.velocity;

            const Vector3R B_p = interpolateB(fieldB, shape, shape05);

            const Vector3R b = 0.5 * dt * q_m * B_p;

            const double betaI = qp * mpw / (1.0 + b.squared());

            const Vector3R I_p =
                betaI * (velocity + velocity.cross(b) + b * velocity.dot(b));

            decompose_current(shape, shape05, I_p, cur_buff);
        }
        auto [start_x, start_y, start_z] = shape.start_.split();
        flush_current_buffer(fieldJ, cur_buff, start_x - 1, start_y - 1,
                             start_z - 1);
    }
}

void calculate_current(const ParticlesArray& particles, Field3d& fieldJ) {
    const double qp = particles.charge;
    const double mpw = particles.mpw();
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
            const Vector3R cell_coord =
                domain.to_cell_coordinates(particle.coord);
            shape.fill_from_normalized(cell_coord);
            shape05.fill_from_normalized(cell_coord, Vector3R(0.5, 0.5, 0.5));
            // TODO: use only current velocity
            const Vector3R velocity =
                0.5 * (particle.velocity + particle.initVelocity);

            const Vector3R I_p = qp * mpw * velocity;

            decompose_current(shape, shape05, I_p, cur_buff);
        }
        auto [start_x, start_y, start_z] = shape.start_.split();
        flush_current_buffer(fieldJ, cur_buff, start_x - 1, start_y - 1,
                             start_z - 1);
    }
}

void predict_current_impl_ngp(const ParticlesArray& particles,
                              const Field3d& fieldB, Field3d& fieldJ,
                              const double dt) {
    const double qp = particles.charge;
    const double mpw = particles.mpw();
    const double q_m = qp / particles.mass();
    const auto& domain = particles.get_domain();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < particles.size(); ++pk) {
        for (auto& particle : particles.particlesData(pk)) {
            const Vector3R cell_coord =
                domain.to_cell_coordinates(particle.coord);

            const Vector3R velocity = particle.velocity;

            double x = cell_coord.x() + GHOST_CELLS;
            double y = cell_coord.y() + GHOST_CELLS;
            double z = cell_coord.z() + GHOST_CELLS;
            double x05 = x - 0.5;
            double y05 = y - 0.5;
            double z05 = z - 0.5;

            const auto ix = ngp(x);
            const auto iy = ngp(y);
            const auto iz = ngp(z);
            const auto ix05 = ngp(x05);
            const auto iy05 = ngp(y05);
            const auto iz05 = ngp(z05);

            Vector3R B_p;
            B_p.x() = fieldB(ix, iy05, iz05, 0);
            B_p.y() = fieldB(ix05, iy, iz05, 1);
            B_p.z() = fieldB(ix05, iy05, iz, 2);

            const Vector3R b = 0.5 * dt * q_m * B_p;

            const double betaI = qp * mpw / (1.0 + b.squared());

            const Vector3R I_p =
                betaI * (velocity + velocity.cross(b) + b * velocity.dot(b));

#pragma omp atomic update
            fieldJ(ix05, iy, iz, 0) += I_p.x();
#pragma omp atomic update
            fieldJ(ix, iy05, iz, 1) += I_p.y();
#pragma omp atomic update
            fieldJ(ix, iy, iz05, 2) += I_p.z();
        }
    }
}

void predict_current(const ParticlesArray& particles, const Field3d& fieldB,
                     Field3d& fieldJ, const double dt, ShapeType type) {
    if (particles.is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            predict_current_impl_ngp(particles, fieldB, fieldJ, dt);
            break;
        case ShapeType::Linear:
            predict_current_impl_linear(particles, fieldB, fieldJ, dt);
            break;
        case ShapeType::Quadratic:
            std::cout
                << "Predict current for quadratic shape not implemented\n";
            exit(-1);
    }
}

Vector3R calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ,
                           const Grid& grid,
                           const BoundaryConditionHandler& bc) {
    Vector3R potE = Vector3R(0,0,0);
  IndexRange range = bc.active_range(grid);

    for (auto i = range.start.x(); i < range.end.x(); ++i) {
        for (auto j = range.start.y(); j < range.end.y(); ++j) {
            for (auto k = range.start.z(); k < range.end.z(); ++k) {
                Vector3R E = Vector3R(fieldE(i, j, k, 0), fieldE(i, j, k, 1),
                                    fieldE(i, j, k, 2));
                Vector3R J = Vector3R(fieldJ(i, j, k, 0), fieldJ(i, j, k, 1),
                                    fieldJ(i, j, k, 2));
                potE += Vector3R(J.x() * E.x(), J.y() * E.y(), J.z() * E.z() );
            }
        }
    }
    return potE;
}

}   // namespace algorithmsECSIM
