#include "algorithms_ecsim.h"

#include "coupling/gather.h"
#include "particles/boris_pusher.h"

namespace algorithmsECSIM {

void predict_velocity_impl_linear(ParticlesArray& particles, const Field3d& fieldEp, const Field3d& fieldB,
                                  const double dt) {
    RECORD_TIMER;
    const double qm = particles.charge() / particles.mass();
    const auto& domain = particles.get_domain();
    require_gather_grid(fieldEp, domain.grid);
    require_gather_grid(fieldB, domain.grid);
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < particles.size(); ++k) {
        for (Particle& particle : particles(k)) {
            const Vector3R normalized_coord = domain.to_cell_coordinates(particle.coord);
            const Stencil3D<2> no = make_point_stencil<ShapeType::Linear>(normalized_coord);
            const Stencil3D<2> sh = make_point_stencil<ShapeType::Linear>(normalized_coord, Vector3R::Constant(0.5));
            const Vector3R E_p = interpolateE(fieldEp, no, sh);
            const Vector3R B_p = interpolateB(fieldB, no, sh);
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

template void predict_velocity<ShapeType::Linear>(ParticlesArray& particles, const Field3d& fieldEp,
                                                  const Field3d& fieldB, const double dt);

Vector3R calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ, const Grid& grid,
                           const FieldBoundaryConditions& bc) {
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
