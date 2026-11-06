#include "coupling/current.h"

#include "coupling/gather.h"

namespace coupling {
void deposit_point_current(const ParticlesArray& particles, Field3d& fieldJ) {
    RECORD_TIMER;

    const double qp = particles.charge();
    const double mpw = particles.mpw();
    const Domain& domain = particles.get_domain();

#pragma omp parallel for schedule(dynamic, particles.omp_granularity())
    for (int pk = 0; pk < particles.size(); ++pk) {
        if (particles(pk).empty()) {
            continue;
        }
        const Vector3R first_coord = domain.to_cell_coordinates(particles(pk).front().coord);
        const Vector3I buffer_origin = make_point_stencil<ShapeType::Linear>(first_coord).origin - Vector3I(1);
        CurrentBuffer<3> cur_buff;
        cur_buff.zero();
        for (const Particle& particle : particles(pk)) {
            const Vector3R cell_coord = domain.to_cell_coordinates(particle.coord);
            const auto shape = make_point_stencil<ShapeType::Linear>(cell_coord);
            const auto shape05 = make_point_stencil<ShapeType::Linear>(cell_coord, Vector3R(0.5));
            const Vector3R velocity = particle.velocity;

            const Vector3R I_p = qp * mpw * velocity;

            decompose_point_current(shape, shape05, I_p, cur_buff, buffer_origin);
        }
        flush_point_current(fieldJ, cur_buff, buffer_origin, domain.grid);
    }
}

}   // namespace coupling
