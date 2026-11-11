#include "drift_kinetic_moments.h"

#include <stdexcept>

#include "coupling/moments.h"
#include "coupling/stencil.h"

namespace driftKinetic {
namespace {

template <ShapeType ST>
void deposit_point_density(const DriftKineticParticles& particles, Field3d& density) {
    const Domain& domain = particles.domain();
#pragma omp parallel for schedule(dynamic, 64)
    for (int cell = 0; cell < particles.size(); ++cell)
        for (const GuidingCenter& p : particles(cell)) {
            const double weight = p.charge * p.weight;   // as in the current deposit
            const auto stencil = make_point_stencil<ST>(domain.to_cell_coordinates(p.position));
            constexpr int N = ST == ShapeType::Linear ? 2 : 3;
            const Vector3I start = to_storage_index(stencil.origin);
            for (int x = 0; x < N; ++x)
                for (int y = 0; y < N; ++y)
                    for (int z = 0; z < N; ++z) {
                        const double value = weight * stencil(x, Axis::X) * stencil(y, Axis::Y) * stencil(z, Axis::Z);
#pragma omp atomic update
                        density(start.x() + x, start.y() + y, start.z() + z, 0) += value;
                    }
        }
}

}   // namespace

void GuidingCenterMoments::deposit_density(Field3d& density, ShapeType shape) const {
    density.setZero();
    if (shape == ShapeType::Linear)
        deposit_point_density<ShapeType::Linear>(particles(), density);
    else if (shape == ShapeType::Quadratic)
        deposit_point_density<ShapeType::Quadratic>(particles(), density);
    else
        throw std::invalid_argument("Guiding-centre density supports Linear and Quadratic shapes");
}

void GuidingCenterMoments::deposit_second_moment(Field3d& moment, int a, int b) const {
    const DriftKineticParticles& sp = particles();
    const Domain& domain = sp.domain();
    const Field3d magnetic = magnetic_ + external_magnetic_;
    const Vector3R center =
        domain.geom.use_cylinder ? domain.geom.cyl_center : 0.5 * (domain.geom.box_min + domain.geom.box_max);
    const Vector3R unit[3] = {Vector3R(1, 0, 0), Vector3R(0, 1, 0), Vector3R(0, 0, 1)};
    moment.setZero();
#pragma omp parallel for schedule(dynamic, 64)
    for (int cell = 0; cell < sp.size(); ++cell)
        for (const GuidingCenter& p : sp(cell)) {
            const Vector3R B = coupling_.magnetic_at(magnetic, p.position);
            const double lenB = B.length();
            // Cylindrical components of b; the basis (r, phi, z) is orthonormal.
            const Vector3R bc = coupling::cylindrical_velocity(p.position, B / lenB, center);
            const double perpendicular = p.mu * lenB / p.mass;   // <v_perp^2> / 2 per perpendicular direction
            const double delta = unit[a].dot(unit[b]);
            const double tensor =
                p.mass * p.weight *
                (p.v_parallel * p.v_parallel * bc[a] * bc[b] + perpendicular * (delta - bc[a] * bc[b]));
            const auto stencil = make_point_stencil<ShapeType::Linear>(domain.to_cell_coordinates(p.position));
            const Vector3I start = to_storage_index(stencil.origin);
            for (int x = 0; x < 2; ++x)
                for (int y = 0; y < 2; ++y)
                    for (int z = 0; z < 2; ++z) {
                        const double value = tensor * stencil(x, Axis::X) * stencil(y, Axis::Y) * stencil(z, Axis::Z);
#pragma omp atomic update
                        moment(start.x() + x, start.y() + y, start.z() + z, 0) += value;
                    }
        }
}

}   // namespace driftKinetic
