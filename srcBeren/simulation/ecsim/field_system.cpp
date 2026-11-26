#include "field_system.h"

#include <algorithm>
#include <cmath>
#include <limits>

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

    domain_ = domain;
    periodic_ = {boundaries.is_periodic(0), boundaries.is_periodic(1), boundaries.is_periodic(2)};
    covered_cells_.resize(domain.size());
    node_cells_.resize(domain.size());
    dilate_scratch_.resize(domain.size());

    cover_margin_ = min_cover_margin;
    // The first cover is not a frequent rebuild and keeps the minimal margin.
    covers_since_rebuild_ = frequent_rebuild_covers;
    // No cell is covered: L gets the pattern of I-M. It is needed before the first
    // cover: without charged particles cover_cells() never rebuilds, and the
    // predictor is I-M.
    build_response(std::vector<int8_t>(unknowns(), 0));
}

namespace {
// result = flags dilated by `radius` cells along one axis, clamped to the storage.
void dilate_axis(const Array3D<int8_t>& flags, int axis, int radius, Array3D<int8_t>& result) {
    const Vector3I size = flags.size();
    const int stride = axis == 0 ? size.y() * size.z() : axis == 1 ? size.z() : 1;
    const int cells = flags.capacity();

#pragma omp parallel for schedule(static)
    for (int cell = 0; cell < cells; ++cell) {
        const int coordinate = (cell / stride) % size[axis];
        const int lo = std::max(0, coordinate - radius) - coordinate;
        const int hi = std::min(size[axis] - 1, coordinate + radius) - coordinate;
        int8_t flagged = 0;
        for (int d = lo; d <= hi && !flagged; ++d) flagged = flags(cell + d * stride);
        result(cell) = flagged;
    }
}

// result = cells within `radius` of a flagged cell along every axis (a box
// dilation), clamped to the storage. One pass per axis.
void dilate(const Array3D<int8_t>& flags, int radius, Array3D<int8_t>& scratch, Array3D<int8_t>& result) {
    RECORD_TIMER;
    dilate_axis(flags, 0, radius, result);
    dilate_axis(result, 1, radius, scratch);
    dilate_axis(scratch, 2, radius, result);
}
}   // namespace

bool EcsimFieldSystem::cover_cells(const Array3D<int8_t>& occupied) {
    RECORD_TIMER;
    assert(occupied.size() == covered_cells_.size() && "one flag per storage cell");
    ++covers_since_rebuild_;
    const int cells = occupied.capacity();
    bool uncovered = false;
#pragma omp parallel for schedule(static) reduction(|| : uncovered)
    for (int cell = 0; cell < cells; ++cell) uncovered = uncovered || (occupied(cell) && !covered_cells_(cell));
    if (!uncovered)
        return false;

    timer::commonTimer timerCover("cover uncovered cells");
    // Margin: doubled after a rebuild that follows the previous one within
    // frequent_rebuild_covers steps (the plasma grows), halved otherwise.
    if (covers_since_rebuild_ <= frequent_rebuild_covers)
        cover_margin_ = std::min(2 * cover_margin_, max_cover_margin);
    else
        cover_margin_ = std::max(min_cover_margin, cover_margin_ / 2);
    covers_since_rebuild_ = 0;

    dilate(occupied, cover_margin_, dilate_scratch_, covered_cells_);
    int64_t covered = 0;
#pragma omp parallel for schedule(static) reduction(+ : covered)
    for (int cell = 0; cell < cells; ++cell) covered += covered_cells_(cell);

    if (covered > full_cover_fraction * cells) {
#pragma omp parallel for schedule(static)
        for (int cell = 0; cell < cells; ++cell) covered_cells_(cell) = 1;
    }

    // A cell writes the rows of nodes from one below to one above it on each
    // axis; the images of a node give rows to its owner, so images of different
    // storage nodes may flag the same row.
    dilate(covered_cells_, 1, dilate_scratch_, node_cells_);
    const Vector3I size = domain_.size();
    std::vector<int8_t> active_rows(unknowns(), 0);

#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int k = 0; k < size.z(); ++k) {
                if (!node_cells_(i, j, k))
                    continue;

                for (int d = 0; d < 3; ++d) {
                    const int unknown = unknown_index[domain_.vind(i, j, k, d)];
                    if (unknown >= 0) {
#pragma omp atomic write
                        active_rows[unknown] = 1;
                    }
                }
            }
    timerCover.finish();

    build_response(active_rows);
    return true;
}

void EcsimFieldSystem::build_response(const std::vector<int8_t>& active_rows) {
    RECORD_TIMER;
    L = ecsim_response::response_pattern(domain_, periodic_, unknown_index, storage_index, active_rows, &IM);
    response_layout = ecsim_response::response_layout(L, domain_, unknown_index, storage_index);

    // Nonzero I-M entries by their positions in L: both rows are sorted, and the
    // columns of an I-M row are a subset of the columns of the L row.
    const int n = unknowns();
    im_start_ = SmartPtr<int>(n + 1);
    im_start_[0] = 0;

#pragma omp parallel for schedule(dynamic, 1024)
    for (int row = 0; row < n; ++row) {
        int nonzero = 0;
        for (Operator::InnerIterator it(IM, row); it; ++it) nonzero += it.value() != 0.0;
        im_start_[row + 1] = nonzero;
    }

    for (int row = 0; row < n; ++row) im_start_[row + 1] += im_start_[row];

    im_position_ = SmartPtr<int>(im_start_[n]);
    im_value_ = SmartPtr<double>(im_start_[n]);
    const VectorView<const int> columns(L.innerIndexPtr(), L.nonZeros());

#pragma omp parallel for schedule(dynamic, 1024)
    for (int row = 0; row < n; ++row) {
        int p = L.outerIndexPtr()[row];
        int k = im_start_[row];
        for (Operator::InnerIterator it(IM, row); it; ++it) {
            while (columns[p] != it.col()) {
                ++p;
                assert(p < L.outerIndexPtr()[row + 1] && "the pattern of L contains I-M");
            }
            if (it.value() != 0.0) {
                im_position_[k] = p;
                im_value_[k] = it.value();
                ++k;
            }
        }
    }
}

void EcsimFieldSystem::zero_response() {
    const int entries = L.nonZeros();
    RECORD_TIMER_PARAMS(int64_t(entries) * int64_t(sizeof(double)), timer::MeasureUnit::byte);
    const VectorView<double> values(L.valuePtr(), entries);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < entries; ++i) values[i] = 0.0;
}

// An entry of I-M+L is kept when I-M has it or its response exceeds
// LMAT_VALUE_TOLERANCE, as the per-step assembly kept before the fixed pattern.
// Rows are counted and filled in chunks taken dynamically: a row with plasma
// response keeps up to 123 entries, a vacuum row only its compact I-M row, so
// equal row ranges per thread are unbalanced when the plasma occupies part of
// the box.
const Operator& EcsimFieldSystem::predictor_matrix() {
    RECORD_TIMER;
    constexpr int chunk_rows = 2048;
    const int n = unknowns();
    assert(L.rows() == n && im_start_[n] <= L.nonZeros() && "L must keep its fixed pattern");

    const VectorView<const int> outer(L.outerIndexPtr(), n + 1);
    const VectorView<const int> columns(L.innerIndexPtr(), L.nonZeros());
    const VectorView<const double> response(L.valuePtr(), L.nonZeros());
    const VectorView<const int> im_start(im_start_.get(), n + 1);
    const VectorView<const int> im_position(im_position_.get(), im_start[n]);
    const VectorView<const double> im_value(im_value_.get(), im_start[n]);

    if (predictor_.rows() != n)
        predictor_.resize(n, n);

    const VectorView<int> starts(predictor_.outerIndexPtr(), n + 1);
    const int chunks = (n + chunk_rows - 1) / chunk_rows;
    if (chunk_offsets_.size != chunks + 1)
        chunk_offsets_ = SmartPtr<int>(chunks + 1);
    const VectorView<int> chunk_offsets(chunk_offsets_.get(), chunks + 1);

    // Kept entries: the response entries above the tolerance and the I-M
    // entries whose response is below it. starts[row + 1] holds the count up to
    // the end of the row within its chunk; the chunk offset is added on filling.
    timer::commonTimer timerCount("count kept entries");
#pragma omp parallel for schedule(dynamic, 1)
    for (int chunk = 0; chunk < chunks; ++chunk) {
        const int begin = chunk * chunk_rows;
        const int end = std::min(n, begin + chunk_rows);
        int kept = 0;

        for (int row = begin; row < end; ++row) {
            for (int p = outer[row]; p < outer[row + 1]; ++p) kept += std::abs(response[p]) > LMAT_VALUE_TOLERANCE;
            for (int k = im_start[row]; k < im_start[row + 1]; ++k)
                kept += std::abs(response[im_position[k]]) <= LMAT_VALUE_TOLERANCE;
            starts[row + 1] = kept;
        }
        chunk_offsets[chunk + 1] = kept;
    }
    timerCount.finish();

    chunk_offsets[0] = 0;
    for (int chunk = 0; chunk < chunks; ++chunk) chunk_offsets[chunk + 1] += chunk_offsets[chunk];
    starts[0] = 0;

    timer::commonTimer timerResize("predictor_.resizeNonZeros(...)");
    predictor_.resizeNonZeros(chunk_offsets[chunks]);
    timerResize.finish();

    // The I-M entries of a row come in the order of its L entries.
    const VectorView<int> kept_columns(predictor_.innerIndexPtr(), chunk_offsets[chunks]);
    const VectorView<double> values(predictor_.valuePtr(), chunk_offsets[chunks]);

    timer::commonTimer timerFill("fill kept entries");
#pragma omp parallel for schedule(dynamic, 1)
    for (int chunk = 0; chunk < chunks; ++chunk) {
        const int begin = chunk * chunk_rows;
        const int end = std::min(n, begin + chunk_rows);
        const int offset = chunk_offsets[chunk];
        int q = offset;
        int k = im_start[begin];

        for (int row = begin; row < end; ++row) {
            starts[row + 1] += offset;
            for (int p = outer[row]; p < outer[row + 1]; ++p) {
                const bool in_im = k < im_start[row + 1] && im_position[k] == p;
                if (in_im || std::abs(response[p]) > LMAT_VALUE_TOLERANCE) {
                    kept_columns[q] = columns[p];
                    values[q] = in_im ? im_value[k] + response[p] : response[p];
                    ++q;
                }
                k += in_im;
            }
            assert(q == starts[row + 1] && "fill matches count");
        }
    }
    timerFill.finish();
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
