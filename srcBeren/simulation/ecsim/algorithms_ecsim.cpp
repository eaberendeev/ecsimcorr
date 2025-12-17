#include "interpolation.h"
#include "simulation_ecsim.h"

namespace algorithmsECSIM {

void predict_velocity_impl_linear(ParticlesArray& particles,
                                  const Field3d& fieldEp, const Field3d& fieldB,
                                  const Domain& domain, const double dt) {
    const double qm = particles.charge / particles.mass();
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < particles.size(); ++k) {
        for (auto& particle : particles.particlesData(k)) {
            const auto normalized_coord =
                particles.to_cell_coordinates(particle.coord);
            const Vector3R E_p = interpolateE_linear(fieldEp, normalized_coord);
            const Vector3R B_p = interpolateB_linear(fieldB, normalized_coord);
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);
        }
    }
}

void predict_velocity_impl_ngp(ParticlesArray& particles,
                               const Field3d& fieldEp, const Field3d& fieldB,
                               [[maybe_unused]] const Domain& domain,
                               const double dt) {
    const double qm = particles.charge / particles.mass();

#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < particles.size(); ++k) {
        for (auto& particle : particles.particlesData(k)) {
            const auto normalized_coord = particles.to_cell_coordinates(particle.coord);
            const Vector3R E_p = interpolateE_ngp(fieldEp, normalized_coord);
            const Vector3R B_p = interpolateB_ngp(fieldB, normalized_coord);
            borisPusher::update_vEB(particle, qm, E_p, B_p, dt);
        }
    }
}

void predict_velocity(ParticlesArray& particles, const Field3d& fieldEp,
                      const Field3d& fieldB, const Domain& domain,
                      const double dt, ShapeType type) {
    if (particles.is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            predict_velocity_impl_ngp(particles, fieldEp, fieldB, domain, dt);
            break;
        case ShapeType::Linear:
            predict_velocity_impl_linear(particles, fieldEp, fieldB, domain,
                                         dt);
            break;
        case ShapeType::Quadratic:
            std::cout
                << "Predict velocity for quadratic shape not implemented\n";
            exit(-1);
    }
}

void predict_current_impl_linear(const ParticlesArray& particles,
                                 const Field3d& fieldB, Field3d& fieldJ,
                                 const Domain& domain, const double dt) {
    constexpr auto SMAX = SHAPE_SIZE;
    const double qp = particles.charge;
    const double mpw = particles.mpw();
    const double q_m = qp / particles.mass();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < particles.size(); ++pk) {
        alignas(64) double wx[SMAX], wy[SMAX], wz[SMAX];
        alignas(64) double wx05[SMAX], wy05[SMAX], wz05[SMAX];
        // TODO: use shape
        // ParticleShape<Shape, 3> shape;
        for (auto& particle : particles.particlesData(pk)) {
            const Vector3R cell_coord = particles.to_cell_coordinates(particle.coord);
           // shape.fill_from_normalized(cell_coord, GHOST_CELLS);
            const Vector3R velocity = particle.velocity;

            double x = cell_coord.x() + GHOST_CELLS;
            double y = cell_coord.y() + GHOST_CELLS;
            double z = cell_coord.z() + GHOST_CELLS;
            double x05 = x - 0.5;
            double y05 = y - 0.5;
            double z05 = z - 0.5;

            const auto ix = int(x);
            const auto iy = int(y);
            const auto iz = int(z);
            const auto ix05 = int(x05);
            const auto iy05 = int(y05);
            const auto iz05 = int(z05);

            wx[1] = (x - ix);
            wx[0] = 1 - wx[1];
            wy[1] = (y - iy);
            wy[0] = 1 - wy[1];
            wz[1] = (z - iz);
            wz[0] = 1 - wz[1];

            wx05[1] = (x05 - ix05);
            wx05[0] = 1 - wx05[1];
            wy05[1] = (y05 - iy05);
            wy05[0] = 1 - wy05[1];
            wz05[1] = (z05 - iz05);
            wz05[0] = 1 - wz05[1];

            const Vector3R B_p = interpolateB_linear(fieldB, cell_coord);

            const Vector3R b = 0.5 * dt * q_m * B_p;

            const double betaI = qp * mpw / (1.0 + b.squared());

            const Vector3R I_p =
                betaI * (velocity + velocity.cross(b) + b * velocity.dot(b));

            for (int nx = 0; nx < SMAX; ++nx) {
                const int i = ix + nx;
                const int i05 = ix05 + nx;
                for (int ny = 0; ny < SMAX; ++ny) {
                    const int j = iy + ny;
                    const int j05 = iy05 + ny;
                    for (int nz = 0; nz < SMAX; ++nz) {
                        const int k = iz + nz;
                        const int k05 = iz05 + nz;
                        double sx = wx05[nx] * wy[ny] * wz[nz];
                        double sy = wx[nx] * wy05[ny] * wz[nz];
                        double sz = wx[nx] * wy[ny] * wz05[nz];
#pragma omp atomic update
                        fieldJ(i05, j, k, 0) += sx * I_p.x();
#pragma omp atomic update
                        fieldJ(i, j05, k, 1) += sy * I_p.y();
#pragma omp atomic update
                        fieldJ(i, j, k05, 2) += sz * I_p.z();
                    }
                }
            }
        }
    }
}

void predict_current_impl_ngp(const ParticlesArray& particles, const Field3d& fieldB,
                              Field3d& fieldJ,
                              [[maybe_unused]] const Domain& domain,
                              const double dt) {
    const double qp = particles.charge;
    const double mpw = particles.mpw();
    const double q_m = qp / particles.mass();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < particles.size(); ++pk) {
        for (auto& particle : particles.particlesData(pk)) {
            const Vector3R cell_coord =
                particles.to_cell_coordinates(particle.coord);

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
                     Field3d& fieldJ, const Domain& domain, const double dt,
                     ShapeType type) {
    if (particles.is_neutral())
        return;

    switch (type) {
        case ShapeType::NGP:
            predict_current_impl_ngp(particles, fieldB, fieldJ, domain, dt);
            break;
        case ShapeType::Linear:
            predict_current_impl_linear(particles, fieldB, fieldJ, domain, dt);
            break;
        case ShapeType::Quadratic:
            std::cout
                << "Predict current for quadratic shape not implemented\n";
            exit(-1);
    }
}

}
