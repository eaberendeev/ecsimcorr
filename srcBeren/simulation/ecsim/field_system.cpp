#include "field_system.h"

#include "runtime/profiler.h"

namespace {
Vector3I owner_node(const Grid& grid, const FieldBoundaryConditions& boundaries, Vector3I node) {
    for (int axis = 0; axis < 3; ++axis)
        if (boundaries.is_periodic(axis))
            node[axis] = grid.periodic_owner(node[axis], axis);
    return node;
}
}   // namespace

void EcsimFieldSystem::init(const MaxwellStencils& stencils, const Domain& domain, double dt,
                            const FieldBoundaryConditions& boundaries) {
    RECORD_TIMER;
    const Grid& grid = domain.grid;
    const int n = 3 * domain.total_size();

    // Owners are enumerated in storage order, so the unknowns keep the node layout.
    unknown_index.assign(n, -1);
    storage_index.clear();
    for (int i = 0; i < grid.size(0); ++i)
        for (int j = 0; j < grid.size(1); ++j)
            for (int k = 0; k < grid.size(2); ++k) {
                if (owner_node(grid, boundaries, Vector3I(i, j, k)) != Vector3I(i, j, k))
                    continue;
                for (int d = 0; d < 3; ++d) {
                    unknown_index[domain.vind(i, j, k, d)] = std::ssize(storage_index);
                    storage_index.push_back(domain.vind(i, j, k, d));
                }
            }
    for (int i = 0; i < grid.size(0); ++i)
        for (int j = 0; j < grid.size(1); ++j)
            for (int k = 0; k < grid.size(2); ++k) {
                const Vector3I owner = owner_node(grid, boundaries, Vector3I(i, j, k));
                for (int d = 0; d < 3; ++d)
                    unknown_index[domain.vind(i, j, k, d)] =
                        unknown_index[domain.vind(owner.x(), owner.y(), owner.z(), d)];
            }

    Operator identity(n, n);
    identity.setIdentity();
    M = -0.25 * dt * dt * stencils.curlB * stencils.curlE;
    Operator storage_IM = identity - M;
    storage_IM.makeCompressed();
    // Row-only pinning: response rows are zeroed separately before A=IM+L.
    // Both predictor and endpoint RHS use the same prescribed electric values.
    impose_dirichlet_rows(storage_IM, boundaries.electric_dirichlet(domain));

    // Owner rows; the columns of all images of a node act on its one unknown.
    std::vector<Trip> trips;
    trips.reserve(storage_IM.nonZeros());
    for (int row = 0; row < unknowns(); ++row)
        for (Operator::InnerIterator it(storage_IM, storage_index[row]); it; ++it)
            trips.emplace_back(row, unknown_index[it.col()], it.value());
    IM.resize(unknowns(), unknowns());
    IM.setFromTriplets(trips.begin(), trips.end());
    IM.makeCompressed();
    L.resize(unknowns(), unknowns());
}

Field3d EcsimFieldSystem::restrict(const Field3d& field) const {
    Field3d result(unknowns());
#pragma omp parallel for schedule(static)
    for (int u = 0; u < unknowns(); ++u) result[u] = field[storage_index[u]];
    return result;
}

void EcsimFieldSystem::prolong(const Field3d& unknowns, Field3d& field) const {
#pragma omp parallel for schedule(static)
    for (int s = 0; s < std::ssize(unknown_index); ++s) field[s] = unknowns[unknown_index[s]];
}
