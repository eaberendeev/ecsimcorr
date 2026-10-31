#pragma once

#include "ParticlesArray.h"
#include "coupling/esirkepov.h"

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
    ::move_and_calc_current<ST>(sp.particlesData, sp.get_domain(), sp.charge, sp.mpw(), dt, fieldJ, policy);
}

template <ShapeType ST = SHAPE>
inline void drift_and_deposit_current_reference(ParticlesArray& sp, double dt, Field3d& fieldJ,
                                                const PathPolicy& policy = PathPolicy::none()) {
    if (sp.is_neutral())
        return;
    ::move_and_calc_current_reference<ST>(sp.particlesData, sp.get_domain(), sp.charge, sp.mpw(), dt, fieldJ, policy);
}

}   // namespace coupling
