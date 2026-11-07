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
        const Vector3I buffer_origin = deposition_origin(make_point_stencil<ShapeType::Linear>(first_coord).origin);
        CurrentBuffer<3> cur_buff;
        cur_buff.zero();
        for (const Particle& particle : particles(pk)) {
            const Vector3R cell_coord = domain.to_cell_coordinates(particle.coord);
            const Stencil3D<2> shape = make_point_stencil<ShapeType::Linear>(cell_coord);
            const Stencil3D<2> shape05 = make_point_stencil<ShapeType::Linear>(cell_coord, Vector3R(0.5));
            const Vector3R velocity = particle.velocity;

            const Vector3R I_p = qp * mpw * velocity;

            decompose_point_current(shape, shape05, I_p, cur_buff, buffer_origin);
        }
        flush_current_buffer(fieldJ, cur_buff, to_storage_index(buffer_origin));
    }
}

}   // namespace coupling
