#pragma once

#include "coupling/esirkepov.h"
#include "coupling/stencil.h"
#include "particles/particle_array.h"

// Accumulate Yee point weights into a buffer with an explicit logical origin.
// The caller groups particles sharing the same unshifted point stencil origin.
template <int N>
inline void decompose_point_current(const Stencil3D<N>& no, const Stencil3D<N>& sh, const Vector3R& value,
                                    CurrentBuffer<N + 1>& buffer, const Vector3I& buffer_origin) {
    for (int x = 0; x < N; ++x) {
        const int ix = no.node(x, Axis::X) - buffer_origin.x();
        const int ix05 = sh.node(x, Axis::X) - buffer_origin.x();
        assert(ix >= 0 && ix < N + 1 && ix05 >= 0 && ix05 < N + 1);
        for (int y = 0; y < N; ++y) {
            const int iy = no.node(y, Axis::Y) - buffer_origin.y();
            const int iy05 = sh.node(y, Axis::Y) - buffer_origin.y();
            assert(iy >= 0 && iy < N + 1 && iy05 >= 0 && iy05 < N + 1);
            for (int z = 0; z < N; ++z) {
                const int iz = no.node(z, Axis::Z) - buffer_origin.z();
                const int iz05 = sh.node(z, Axis::Z) - buffer_origin.z();
                assert(iz >= 0 && iz < N + 1 && iz05 >= 0 && iz05 < N + 1);
                buffer(ix05, iy, iz, 0) += value.x() * sh(x, Axis::X) * no(y, Axis::Y) * no(z, Axis::Z);
                buffer(ix, iy05, iz, 1) += value.y() * no(x, Axis::X) * sh(y, Axis::Y) * no(z, Axis::Z);
                buffer(ix, iy, iz05, 2) += value.z() * no(x, Axis::X) * no(y, Axis::Y) * sh(z, Axis::Z);
            }
        }
    }
}

// Point-current buffers can have zero padding beyond the upper wall. Keep the
// existing wall clipping here explicitly, independent of logical-node access.
template <int N>
inline void flush_point_current(Field3d& field, const CurrentBuffer<N>& buffer, const Vector3I& logical_origin,
                                const Grid& grid) {
    flush_current_buffer(field, buffer, logical_origin, grid);
}

namespace coupling {
// Accumulate the instantaneous linear Yee current sum(q*mpw*v*S).
// fieldJ must be cleared by the caller before accumulating species.
void deposit_point_current(const ParticlesArray& particles, Field3d& fieldJ);

// Drift coordinates by dt and accumulate the trajectory-averaged current.
// Every segment (including both reflected legs) uses this same dt for its
// normalization. Wall fates are recorded by the path policy and applied by the
// caller before cell redistribution invalidates particle addresses.
template <ShapeType ST = SHAPE>
inline void drift_and_deposit_current(ParticlesArray& sp, double dt, Field3d& fieldJ,
                                      const PathPolicy& policy = PathPolicy::none()) {
    if (sp.is_neutral())
        return;
    ::move_and_calc_current<ST>(sp.cells(), sp.get_domain(), sp.charge(), sp.mpw(), dt, fieldJ, policy);
}

template <ShapeType ST = SHAPE>
inline void drift_and_deposit_current_reference(ParticlesArray& sp, double dt, Field3d& fieldJ,
                                                const PathPolicy& policy = PathPolicy::none()) {
    if (sp.is_neutral())
        return;
    ::move_and_calc_current_reference<ST>(sp.cells(), sp.get_domain(), sp.charge(), sp.mpw(), dt, fieldJ, policy);
}

}   // namespace coupling
