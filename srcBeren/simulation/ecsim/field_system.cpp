#include "field_system.h"

void EcsimFieldSystem::init(const MaxwellStencils& stencils, const Domain& domain, double dt,
                            BoundaryConditionHandler& boundaries) {
    const int n = 3 * domain.total_size();
    L.resize(n, n);
    L.setZero();
    Operator identity(n, n);
    identity.setIdentity();
    M = -0.25 * dt * dt * stencils.curlB * stencils.curlE;
    IM = identity - M;
    IM.makeCompressed();
    // Row-only pinning: response rows are zeroed separately before A=IM+L.
    // Both predictor and endpoint RHS use the same prescribed electric values.
    impose_dirichlet_rows(IM, boundaries.electric_dirichlet(domain));
}
