#include "coupling/current.h"

#include "coupling/gather.h"

namespace coupling {
void deposit_point_current(const ParticlesArray& particles, Field3d& fieldJ) {
    RECORD_TIMER;

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
            const Vector3R cell_coord = domain.to_cell_coordinates(particle.coord);
            shape.fill_from_normalized(cell_coord);
            shape05.fill_from_normalized(cell_coord, Vector3R(0.5, 0.5, 0.5));
            const Vector3R velocity = particle.velocity;

            const Vector3R I_p = qp * mpw * velocity;

            decompose_current(shape, shape05, I_p, cur_buff);
        }
        auto [start_x, start_y, start_z] = shape.start_.split();
        flush_current_buffer(fieldJ, cur_buff, start_x - 1, start_y - 1, start_z - 1);
    }
}

}   // namespace coupling
