#include "field_system.h"

#include <omp.h>

#include <algorithm>
#include <cmath>

#include "include/config.h"
#include "runtime/profiler.h"
#include "simulation/ecsim/response/response.h"
#include "utils/memory.h"

namespace {
Vector3I owner_node(const Grid& grid, const FieldBoundaryConditions& boundaries, Vector3I node) {
    for (int axis = 0; axis < 3; ++axis)
        if (boundaries.is_periodic(axis))
            node[axis] = grid.periodic_node_owner(node[axis], axis);
    return node;
}

bool keeps_entry(double im, double response) {
    return im != 0.0 || std::abs(response) > LMAT_VALUE_TOLERANCE;
}

bool is_identity_row(const Operator& matrix, int row) {
    for (Operator::InnerIterator it(matrix, row); it; ++it)
        if (it.value() != (it.col() == row ? 1.0 : 0.0))
            return false;
    return true;
}
}   // namespace

void EcsimFieldSystem::init(const MaxwellStencils& stencils, const Domain& domain, double dt,
                            const FieldBoundaryConditions& boundaries) {
    RECORD_TIMER;
    const Grid& grid = domain.grid;
    const int n = 3 * domain.total_size();

    Operator identity(n, n);
    identity.setIdentity();
    M = -0.25 * dt * dt * stencils.curlB * stencils.curlE;
    Operator storage_IM = identity - M;
    storage_IM.makeCompressed();
    // Row-only pinning: response rows are zeroed separately before A=IM+L.
    // Both predictor and endpoint RHS use the same prescribed electric values.
    const std::vector<ElectricDirichletEntry>& dirichlet = boundaries.electric_dirichlet(domain);
    impose_dirichlet_rows(storage_IM, dirichlet);
    std::vector<char> prescribed(n, 0);
    for (const ElectricDirichletEntry& entry : dirichlet) prescribed[entry.row] = 1;

    // Owners are enumerated in storage order, so the unknowns keep the node layout.
    std::vector<int> owner(n);
    unknown_index.assign(n, -1);
    storage_index.clear();
    for (int i = 0; i < grid.size(0); ++i)
        for (int j = 0; j < grid.size(1); ++j)
            for (int k = 0; k < grid.size(2); ++k) {
                const Vector3I node = owner_node(grid, boundaries, Vector3I(i, j, k));
                for (int d = 0; d < 3; ++d) {
                    const int s = domain.vind(i, j, k, d);
                    owner[s] = domain.vind(node.x(), node.y(), node.z(), d);
                    if (owner[s] == s && (prescribed[s] || !is_identity_row(storage_IM, s))) {
                        unknown_index[s] = std::ssize(storage_index);
                        storage_index.push_back(s);
                    }
                }
            }
    // A separate pass: an image in the lower ghost layer precedes its owner in storage.
    for (int s = 0; s < n; ++s) unknown_index[s] = unknown_index[owner[s]];

    // Owner rows; the columns of all images of a node act on its one unknown.
    std::vector<Trip> trips, fixed;
    trips.reserve(storage_IM.nonZeros());
    for (int row = 0; row < unknowns(); ++row)
        for (Operator::InnerIterator it(storage_IM, storage_index[row]); it; ++it) {
            if (unknown_index[it.col()] >= 0)
                trips.emplace_back(row, unknown_index[it.col()], it.value());
            else
                fixed.emplace_back(row, it.col(), it.value());
        }
    IM.resize(unknowns(), unknowns());
    IM.setFromTriplets(trips.begin(), trips.end());
    IM.makeCompressed();
    fixed_columns_.resize(unknowns(), n);
    fixed_columns_.setFromTriplets(fixed.begin(), fixed.end());
    fixed_columns_.makeCompressed();

    // L and I-M+L share one pattern: the couplings of a particle united with I-M.
    const std::array<bool, 3> periodic{boundaries.is_periodic(0), boundaries.is_periodic(1), boundaries.is_periodic(2)};
    L = ecsim_response::response_pattern(domain, periodic, unknown_index, storage_index) + IM;
    L.makeCompressed();
    im_on_pattern_.assign(L.valuePtr(), L.valuePtr() + L.nonZeros());
    clear_response();
}

void EcsimFieldSystem::clear_response() {
    const VectorView<double> values(L.valuePtr(), L.nonZeros());
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L.nonZeros(); ++i) values[i] = 0.0;
}

// I-M+L without the pattern entries that no particle reached: an entry stays when
// I-M has it or its response exceeds LMAT_VALUE_TOLERANCE, as the per-step
// assembly kept before the fixed pattern. Each thread counts its rows, then fills
// them from its offset, so the row starts that balance the partitioned SpMV
// follow the kept entries.
const Operator& EcsimFieldSystem::predictor_matrix() {
    RECORD_TIMER;
    assert(std::ssize(im_on_pattern_) == L.nonZeros() && "L must keep its fixed pattern");
    const int n = unknowns();
    const VectorView<const int> outer(L.outerIndexPtr(), n + 1);
    const VectorView<const int> columns(L.innerIndexPtr(), L.nonZeros());
    const VectorView<const double> response(L.valuePtr(), L.nonZeros());
    const VectorView<const double> im(im_on_pattern_.data(), L.nonZeros());
    predictor_.resize(n, n);
    const VectorView<int> starts(predictor_.outerIndexPtr(), n + 1);
    SmartPtr<int> thread_offsets(omp_get_max_threads() + 1);
#pragma omp parallel
    {
        const int threads = omp_get_num_threads();
        const int thread = omp_get_thread_num();
        const int begin = static_cast<int64_t>(n) * thread / threads;
        const int end = static_cast<int64_t>(n) * (thread + 1) / threads;
        int kept = 0;
        for (int row = begin; row < end; ++row) {
            for (int p = outer[row]; p < outer[row + 1]; ++p) kept += keeps_entry(im[p], response[p]);
            starts[row + 1] = kept;
        }
        thread_offsets[thread + 1] = kept;
#pragma omp barrier
#pragma omp single
        {
            thread_offsets[0] = 0;
            for (int t = 0; t < threads; ++t) thread_offsets[t + 1] += thread_offsets[t];
            predictor_.resizeNonZeros(thread_offsets[threads]);
        }
        const VectorView<int> kept_columns(predictor_.innerIndexPtr(), thread_offsets[threads]);
        const VectorView<double> values(predictor_.valuePtr(), thread_offsets[threads]);
        const int offset = thread_offsets[thread];
        int q = offset;
        for (int row = begin; row < end; ++row) {
            starts[row + 1] += offset;
            for (int p = outer[row]; p < outer[row + 1]; ++p)
                if (keeps_entry(im[p], response[p])) {
                    kept_columns[q] = columns[p];
                    values[q] = im[p] + response[p];
                    ++q;
                }
        }
    }
    return predictor_;
}

Field3d EcsimFieldSystem::restrict(const Field3d& field) const {
    RECORD_TIMER;
    Field3d result(unknowns());
#pragma omp parallel for schedule(static)
    for (int u = 0; u < unknowns(); ++u) result[u] = field[storage_index[u]];
    return result;
}

Field3d EcsimFieldSystem::restrict_rhs(const Field3d& rhs) const {
    RECORD_TIMER;
    Field3d result = restrict(rhs);
    if (fixed_columns_.nonZeros() > 0)
        result.data() -= fixed_columns_ * rhs.data();
    return result;
}

void EcsimFieldSystem::prolong(const Field3d& unknowns, const Field3d& fixed, Field3d& field) const {
    RECORD_TIMER;
#pragma omp parallel for schedule(static)
    for (int s = 0; s < std::ssize(unknown_index); ++s)
        field[s] = unknown_index[s] >= 0 ? unknowns[unknown_index[s]] : fixed[s];
}

void EcsimFieldSystem::prolong(const Field3d& unknowns, Field3d& field) const {
    RECORD_TIMER;
#pragma omp parallel for schedule(static)
    for (int s = 0; s < std::ssize(unknown_index); ++s)
        field[s] = unknown_index[s] >= 0 ? unknowns[unknown_index[s]] : 0.0;
}
