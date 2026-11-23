#include "response.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <iterator>

#include "coupling/stencil.h"
#include "particles/particle_array.h"
#include "response_kernels.h"
#include "runtime/profiler.h"
#include "utils/memory.h"

namespace ecsim_response {
namespace {
// Node offsets from a row of component c1 to the columns of component c2 that one
// Linear particle couples: [-1,1] per axis, [-1,2] along c1 and [-2,1] along c2
// when the components differ.
void coupling_box(int c1, int c2, Vector3I& lo, Vector3I& hi) {
    for (int a = 0; a < 3; ++a) {
        lo[a] = (a == c2 && a != c1) ? -2 : -1;
        hi[a] = (a == c1 && a != c2) ? 2 : 1;
    }
}

constexpr int box_extent = 5;   // offsets -2..2 per axis

int box_index(int c2, const Vector3I& offset) {
    return ((c2 * box_extent + offset.x() + 2) * box_extent + offset.y() + 2) * box_extent + offset.z() + 2;
}

constexpr int box_entries = 3 * box_extent * box_extent * box_extent;

// Couplings of a row of component c1 as (dx, dy, dz, c2), in the order of a regular
// row: unknowns follow storage order, which sorts nodes by x, y, z and then the
// component.
std::vector<std::array<int, 4>> ordered_couplings(int c1) {
    std::vector<std::array<int, 4>> keys;
    for (int c2 = 0; c2 < 3; ++c2) {
        Vector3I lo, hi;
        coupling_box(c1, c2, lo, hi);
        for (int dx = lo.x(); dx <= hi.x(); ++dx)
            for (int dy = lo.y(); dy <= hi.y(); ++dy)
                for (int dz = lo.z(); dz <= hi.z(); ++dz) keys.push_back({dx, dy, dz, c2});
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

// Positions of the couplings of a regular row of component c1, indexed by box_index().
std::array<int, box_entries> canonical_positions(int c1) {
    const std::vector<std::array<int, 4>> keys = ordered_couplings(c1);
    std::array<int, box_entries> positions;
    positions.fill(-1);
    for (int p = 0; p < std::ssize(keys); ++p)
        positions[box_index(keys[p][3], Vector3I(keys[p][0], keys[p][1], keys[p][2]))] = p;
    return positions;
}

template <typename RowIdx, typename ColIdx>
void fill_block_positions(const std::array<int, box_entries>& canonical, std::array<int, layout_entries>& position) {
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1)
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1)
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1)
                for (int x2 = 0; x2 < ColIdx::size_x; ++x2)
                    for (int y2 = 0; y2 < ColIdx::size_y; ++y2)
                        for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                            const Vector3I offset(x2 + ColIdx::offset_x - x1 - RowIdx::offset_x,
                                                  y2 + ColIdx::offset_y - y1 - RowIdx::offset_y,
                                                  z2 + ColIdx::offset_z - z1 - RowIdx::offset_z);
                            bool inside = true;
                            for (int a = 0; a < 3; ++a) inside = inside && offset[a] >= -2 && offset[a] <= 2;
                            const int r = RowIdx::calculate(x1, y1, z1);
                            const int c = ColIdx::calculate(x2, y2, z2);
                            position[((RowIdx::dir * block_nodes + r) * 3 + ColIdx::dir) * block_nodes + c] =
                                inside ? canonical[box_index(ColIdx::dir, offset)] : -1;
                        }
}

template <typename RowIdx>
void fill_row_positions(const std::array<int, box_entries>& canonical, std::array<int, layout_entries>& position) {
    fill_block_positions<RowIdx, XIndexer>(canonical, position);
    fill_block_positions<RowIdx, YIndexer>(canonical, position);
    fill_block_positions<RowIdx, ZIndexer>(canonical, position);
}
}   // namespace

ResponseLayout response_layout(const Operator& response, const Domain& domain, const std::vector<int>& unknown_index,
                               const std::vector<int>& storage_index) {
    RECORD_TIMER;
    const std::array<std::array<int, box_entries>, 3> canonical{canonical_positions(0), canonical_positions(1),
                                                                canonical_positions(2)};
    ResponseLayout layout;
    fill_row_positions<XIndexer>(canonical[0], layout.position);
    fill_row_positions<YIndexer>(canonical[1], layout.position);
    fill_row_positions<ZIndexer>(canonical[2], layout.position);

    const Grid& grid = domain.grid;
    const int n = response.rows();
    int full_row = 0;
    for (int c2 = 0; c2 < 3; ++c2) {
        Vector3I lo, hi;
        coupling_box(0, c2, lo, hi);
        full_row += (hi.x() - lo.x() + 1) * (hi.y() - lo.y() + 1) * (hi.z() - lo.z() + 1);
    }

    // A row is irregular when its columns differ from the 123 of an interior node
    // in storage order:
    // - the node is within two nodes of a wall or of the end of the field domain:
    //   neighbors outside the storage or fixed components drop out, the row is shorter;
    // - the node is within two nodes of a periodic join: neighbors in the ghost layer
    //   are images and take the unknown of their owner at the other end of the box,
    //   so the columns leave storage order (and join on a short periodic axis);
    // - the node is away from the covered cells: the row keeps only its I-M columns
    //   and gets no response.
    layout.regular.assign(n, 0);
#pragma omp parallel for schedule(dynamic, 1024)
    for (int row = 0; row < n; ++row) {
        const int begin = response.outerIndexPtr()[row];
        if (response.outerIndexPtr()[row + 1] - begin != full_row)
            continue;

        const std::array<int, 4> node = grid.pos_vind_range(storage_index[row]);
        const int c1 = node[3];
        bool regular = true;
        for (int c2 = 0; c2 < 3 && regular; ++c2) {
            Vector3I lo, hi;
            coupling_box(c1, c2, lo, hi);
            for (int dx = lo.x(); dx <= hi.x() && regular; ++dx)
                for (int dy = lo.y(); dy <= hi.y() && regular; ++dy)
                    for (int dz = lo.z(); dz <= hi.z() && regular; ++dz) {
                        const Vector3I neighbor(node[0] + dx, node[1] + dy, node[2] + dz);
                        bool stored = true;
                        for (int a = 0; a < 3; ++a) stored = stored && neighbor[a] >= 0 && neighbor[a] < grid.size(a);
                        const int column =
                            stored ? unknown_index[domain.vind(neighbor.x(), neighbor.y(), neighbor.z(), c2)] : -1;
                        const int p = canonical[c1][box_index(c2, Vector3I(dx, dy, dz))];
                        regular = column >= 0 && response.innerIndexPtr()[begin + p] == column;
                    }
        }

        layout.regular[row] = regular;
    }
    return layout;
}

namespace {
// Sorted columns of one row of response_pattern(). The couplings come in the order
// of a regular row, so only rows with periodic joins are sorted; the I-M row, also
// sorted, is merged in.
void pattern_row(int row, const Domain& domain, const std::array<bool, 3>& periodic,
                 const std::vector<int>& unknown_index, const std::vector<int>& storage_index,
                 const std::vector<int8_t>& active_rows, const Operator* united_with,
                 const std::array<std::vector<std::array<int, 4>>, 3>& couplings, std::vector<int>& columns,
                 std::vector<int>& merged) {
    const Grid& grid = domain.grid;
    columns.clear();
    if (active_rows.empty() || active_rows[row]) {
        const std::array<int, 4> node = grid.pos_vind_range(storage_index[row]);
        bool sorted = true;
        for (const std::array<int, 4>& coupling : couplings[node[3]]) {
            Vector3I neighbor(node[0] + coupling[0], node[1] + coupling[1], node[2] + coupling[2]);
            bool stored = true;
            for (int a = 0; a < 3; ++a) {
                if (periodic[a])
                    neighbor[a] = grid.periodic_node_owner(neighbor[a], a);
                stored = stored && neighbor[a] >= 0 && neighbor[a] < grid.size(a);
            }
            const int column =
                stored ? unknown_index[domain.vind(neighbor.x(), neighbor.y(), neighbor.z(), coupling[3])] : -1;
            if (column >= 0) {
                sorted = sorted && (columns.empty() || columns.back() < column);
                columns.push_back(column);
            }
        }
        if (!sorted) {
            std::sort(columns.begin(), columns.end());
            columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
        }
    }
    if (united_with != nullptr) {
        const int begin = united_with->outerIndexPtr()[row];
        const int end = united_with->outerIndexPtr()[row + 1];
        merged.clear();
        std::set_union(columns.begin(), columns.end(), united_with->innerIndexPtr() + begin,
                       united_with->innerIndexPtr() + end, std::back_inserter(merged));
        columns.swap(merged);
    }
}
}   // namespace

Operator response_pattern(const Domain& domain, const std::array<bool, 3>& periodic,
                          const std::vector<int>& unknown_index, const std::vector<int>& storage_index,
                          const std::vector<int8_t>& active_rows, const Operator* united_with) {
    RECORD_TIMER;
    constexpr int chunk_rows = 1024;
    const int n = std::ssize(storage_index);
    assert((active_rows.empty() || std::ssize(active_rows) == n) && "one flag per unknown");
    const std::array<std::vector<std::array<int, 4>>, 3> couplings{ordered_couplings(0), ordered_couplings(1),
                                                                   ordered_couplings(2)};

    Operator pattern(n, n);
    const VectorView<int> starts(pattern.outerIndexPtr(), n + 1);
    // Rows are generated once into chunk buffers, then copied behind their row starts.
    const int chunks = (n + chunk_rows - 1) / chunk_rows;
    std::vector<std::vector<int>> chunk_columns(chunks);
#pragma omp parallel
    {
        std::vector<int> columns, merged;
#pragma omp for schedule(dynamic, 1)
        for (int chunk = 0; chunk < chunks; ++chunk) {
            for (int row = chunk * chunk_rows; row < std::min(n, (chunk + 1) * chunk_rows); ++row) {
                pattern_row(row, domain, periodic, unknown_index, storage_index, active_rows, united_with, couplings,
                            columns, merged);
                starts[row + 1] = std::ssize(columns);
                chunk_columns[chunk].insert(chunk_columns[chunk].end(), columns.begin(), columns.end());
            }
        }
    }

    starts[0] = 0;
    for (int row = 0; row < n; ++row) starts[row + 1] += starts[row];
    pattern.resizeNonZeros(starts[n]);
    const VectorView<int> inner(pattern.innerIndexPtr(), starts[n]);
    const VectorView<double> values(pattern.valuePtr(), starts[n]);
#pragma omp parallel for schedule(dynamic, 1)
    for (int chunk = 0; chunk < chunks; ++chunk) {
        const int offset = starts[chunk * chunk_rows];
        for (int e = 0; e < std::ssize(chunk_columns[chunk]); ++e) {
            inner[offset + e] = chunk_columns[chunk][e];
            values[offset + e] = 0.0;
        }
        std::vector<int>().swap(chunk_columns[chunk]);
    }
    return pattern;
}

void mark_occupied_cells(const ParticlesArray& species, Array3D<int8_t>& occupied) {
    RECORD_TIMER;
    assert(occupied.size() == species.cells().size() && "one flag per storage cell");
#pragma omp parallel for schedule(static)
    for (int cell = 0; cell < species.size(); ++cell)
        if (!species(cell).empty())
            occupied(cell) = 1;
}

namespace {
// Add a sorted row of contributions to the same row of `response`. Cells of one
// color never share a row (accumulate_species), so the update needs no atomics.
void add_row(const RowBlock<3 * block_nodes>& block, Operator& response) {
    const int begin = response.outerIndexPtr()[block.row];
    const int length = response.outerIndexPtr()[block.row + 1] - begin;
    const VectorView<const int> columns(response.innerIndexPtr() + begin, length);
    const VectorView<double> values(response.valuePtr() + begin, length);

    int p = 0;
    for (int e = 0; e < block.nnz; ++e) {
        while (columns[p] < block.columns[e]) ++p;
        assert(columns[p] == block.columns[e] && "response entry outside the fixed pattern");
        values[p] += block.values[e];
    }
}

// Columns are field unknowns: periodic images of a node map to the unknown of its owner,
// which can break the column order and join columns on short periodic axes. Fixed
// components (outside the field domain, index -1) carry no response.
template <typename ColIdx, int DIR, typename Block_t>
void collect_block_row(int rowIdx, int i_cell, int j_cell, int k_cell, const Block_t& block, const Domain& domain,
                       const std::vector<int>& unknown_index, double tolerance, RowBlock<block_nodes>& rowBlock) {
    bool sorted = true;
    for (int x2 = 0; x2 < ColIdx::size_x; ++x2) {
        for (int y2 = 0; y2 < ColIdx::size_y; ++y2) {
            for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                const int colIdx = ColIdx::calculate(x2, y2, z2);
                const double val = block(rowIdx, colIdx, DIR);

                if (std::abs(val) > tolerance) {
                    const int col =
                        unknown_index[domain.vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                                  k_cell + z2 + ColIdx::offset_z, ColIdx::dir)];
                    if (col < 0)
                        continue;
                    sorted = sorted && (rowBlock.nnz == 0 || rowBlock.columns[rowBlock.nnz - 1] < col);
                    rowBlock.push_back_value(col, val);
                }
            }
        }
    }
    if (!sorted)
        rowBlock.sort_and_merge();
}

// Add the couplings of one block row to a regular row at their known positions.
template <int baseDir, typename Block_t>
void add_regular_row(int row, int rowIdx, const Block_t& block, const ResponseLayout& layout, double tolerance,
                     Operator& response) {
    double* values = response.valuePtr() + response.outerIndexPtr()[row];
    const int* position = &layout.position[(baseDir / 3 * block_nodes + rowIdx) * 3 * block_nodes];
    for (int c2 = 0; c2 < 3; ++c2)
        for (int c = 0; c < block_nodes; ++c) {
            const double val = block(rowIdx, c, baseDir + c2);
            const int p = position[c2 * block_nodes + c];
            assert((p >= 0 || val == 0.0) && "a particle couples outside the response pattern");
            if (p >= 0 && std::abs(val) > tolerance)
                values[p] += val;
        }
}

template <typename RowIdx, int baseDir, typename Block_t>
void add_block_rows(int i_cell, int j_cell, int k_cell, const Block_t& block, const Domain& domain,
                    const std::vector<int>& unknown_index, const ResponseLayout& layout, double tolerance,
                    Operator& response) {
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row =
                    unknown_index[domain.vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                              k_cell + z1 + RowIdx::offset_z, RowIdx::dir)];

                if (row < 0)
                    continue;

                const int rowIdx = RowIdx::calculate(x1, y1, z1);
                if (layout.regular[row]) {
                    add_regular_row<baseDir>(row, rowIdx, block, layout, tolerance, response);
                    continue;
                }

                // Rows near walls and periodic joins: collect and merge the
                // couplings, then search their columns.
                RowBlock<block_nodes> rowBuffers[3]{row, row, row};
                collect_block_row<XIndexer, baseDir + 0>(rowIdx, i_cell, j_cell, k_cell, block, domain, unknown_index,
                                                         tolerance, rowBuffers[0]);
                collect_block_row<YIndexer, baseDir + 1>(rowIdx, i_cell, j_cell, k_cell, block, domain, unknown_index,
                                                         tolerance, rowBuffers[1]);
                collect_block_row<ZIndexer, baseDir + 2>(rowIdx, i_cell, j_cell, k_cell, block, domain, unknown_index,
                                                         tolerance, rowBuffers[2]);

                const bool isEmpty = rowBuffers[0].nnz == 0 && rowBuffers[1].nnz == 0 && rowBuffers[2].nnz == 0;
                if (!isEmpty) {
                    RowBlock<3 * block_nodes> tmp;
                    tmp.mergeFromOthers<3>(rowBuffers);
                    add_row(tmp, response);
                }
            }
        }
    }
}

// Response of the particles of one species in one cell, added to `tmpBlock`.
void accumulate_cell_particles(const std::vector<Particle>& currVec, const ParticlesArray& species,
                               const Field3d& fieldB, const Domain& domain, const double dt,
                               [[maybe_unused]] const Vector3I& storage_cell, Block& tmpBlock,
                               PredictedCurrent* current) {
    const Vector3R first_coord = domain.to_cell_coordinates(currVec.front().coord);
    const Vector3I cell_origin = checked_floor_index(first_coord);

    if (currVec.size() >= 64) {
        constexpr int buffSize = 16;
        std::array<std::array<const Particle*, buffSize>, 8> particleBuffers;
        std::array<int, 8> accumulatedCount{};

        for (const Particle& particle : currVec) {
            const Vector3R& coord = particle.coord;
            const Vector3R coordLoc = domain.to_cell_coordinates(coord);
            const Vector3I shifted_origin = linear_half_shifted_origin(coordLoc, cell_origin);
            const int kindX = cell_origin.x() - shifted_origin.x();
            const int kindY = cell_origin.y() - shifted_origin.y();
            const int kindZ = cell_origin.z() - shifted_origin.z();

            const int currKind = (kindX * 2 + kindY) * 2 + kindZ;
            assert(currKind >= 0 && currKind < 8);

            particleBuffers[currKind][accumulatedCount[currKind]] = &particle;
            accumulatedCount[currKind] += 1;

            if (accumulatedCount[currKind] == buffSize) {
                update_Lmat2<buffSize>(&particleBuffers[currKind][0], domain, species.charge(), species.mass(),
                                       species.mpw(), fieldB, dt, tmpBlock, current);
                accumulatedCount[currKind] = 0;
            }
        }

        for (int kind = 0; kind < 8; ++kind) {
            const int accumulated = accumulatedCount[kind];
            int offset = 0;
            if (accumulated >= 8) {
                update_Lmat2<8>(&particleBuffers[kind][offset], domain, species.charge(), species.mass(), species.mpw(),
                                fieldB, dt, tmpBlock, current);
                offset += 8;
            }
            if (accumulated - offset >= 4) {
                update_Lmat2<4>(&particleBuffers[kind][offset], domain, species.charge(), species.mass(), species.mpw(),
                                fieldB, dt, tmpBlock, current);
                offset += 4;
            }
            if (accumulated - offset >= 2) {
                update_Lmat2<2>(&particleBuffers[kind][offset], domain, species.charge(), species.mass(), species.mpw(),
                                fieldB, dt, tmpBlock, current);
                offset += 2;
            }
            if (accumulated - offset >= 1) {
                update_Lmat2(*particleBuffers[kind][offset], domain, species.charge(), species.mass(), species.mpw(),
                             fieldB, dt, tmpBlock, current);
                offset += 1;
            }
            assert(accumulated == offset);
        }
    } else {
        for (const Particle& particle : currVec)
            update_Lmat2(particle, domain, species.charge(), species.mass(), species.mpw(), fieldB, dt, tmpBlock,
                         current);
    }

    assert(to_storage_index(cell_origin) == storage_cell && "particles are grouped by their storage cell");
}

// Add the predicted current of a cell to `field`, the buffer starting one node
// below the cell. Cells of one color never share a node, so no atomics.
void add_predicted_current(const PredictedCurrent& buffer, const Vector3I& storage_cell, Field3d& field) {
    const Vector3I start = storage_cell - Vector3I::Constant(1);
    const Vector3I sizes = field.sizes();
    for (int n = 0; n < 3; ++n) {
        const int ix = start.x() + n;
        for (int m = 0; m < 3; ++m) {
            const int iy = start.y() + m;
            for (int l = 0; l < 3; ++l) {
                const int iz = start.z() + l;
                if (ix < 0 || iy < 0 || iz < 0 || ix >= sizes.x() || iy >= sizes.y() || iz >= sizes.z())
                    continue;

                for (int d = 0; d < 3; ++d) field(ix, iy, iz, d) += buffer(n, m, l, d);
            }
        }
    }
}

// Response of all charged species of one cell, summed into one block and added
// to the rows of that block once. Cells without particles add nothing. With
// `current`, the predicted current of the cell is added to it.
void accumulate_cell(const std::vector<const ParticlesArray*>& charged, const Field3d& fieldB, const Domain& domain,
                     const double dt, const std::vector<int>& unknown_index, const ResponseLayout& layout,
                     const Vector3I& storage_cell, Block& tmpBlock, PredictedCurrent& cellCurrent, Operator& response,
                     Field3d* current) {
    const auto [i, j, k] = storage_cell.split();
    bool isEmpty = true;
    for (const ParticlesArray* species : charged) {
        const std::vector<Particle>& currVec = (*species)(i, j, k);
        if (currVec.empty())
            continue;

        if (isEmpty) {
            tmpBlock.setZero();
            cellCurrent.zero();
        }
        isEmpty = false;
        accumulate_cell_particles(currVec, *species, fieldB, domain, dt, storage_cell, tmpBlock,
                                  current ? &cellCurrent : nullptr);
    }
    if (isEmpty)
        return;

    add_block_rows<XIndexer, 0>(i, j, k, tmpBlock, domain, unknown_index, layout, response_value_tolerance, response);
    add_block_rows<YIndexer, 3>(i, j, k, tmpBlock, domain, unknown_index, layout, response_value_tolerance, response);
    add_block_rows<ZIndexer, 6>(i, j, k, tmpBlock, domain, unknown_index, layout, response_value_tolerance, response);
    if (current)
        add_predicted_current(cellCurrent, storage_cell, *current);
}

// Tiles along one axis: `width` interior cells each (at least two), the rest
// joined to the last tile, ghost cells to the end tiles. Returns the tile
// starts in storage cells and the end of the storage range.
std::vector<int> axis_tiles(int storage_cells, int interior_cells, int width = 2) {
    assert(width >= 2 && "same-color tiles must be two cells apart");
    const int tiles = std::max(1, interior_cells / width);
    std::vector<int> starts(tiles + 1);
    starts[0] = 0;
    for (int t = 1; t < tiles; ++t) starts[t] = GHOST_CELLS + width * t;
    starts[tiles] = storage_cells;
    return starts;
}

// Tiles alternate two colors. A cell writes the rows of nodes up to one cell
// away, and across a periodic boundary the last tile is next to the first one
// (color 0). With an even tile count the last tile has color 1; with an odd
// count it would have color 0 like the first, so it gets a third color.
int tile_color(int tile, int tiles) {
    const bool odd_last = tiles % 2 == 1 && tiles > 1 && tile == tiles - 1;
    return odd_last ? 2 : tile % 2;
}

// Tiles of one color per thread wanted for load balance: in a cylinder or a
// sparse plasma a large share of tiles is empty.
constexpr int min_tiles_per_thread = 8;
}   // namespace

// Cells are taken in tiles of two cells along x and y; along z a tile takes
// the full column, or a part of it when x and y give too few tiles of one color
// for the threads (min_tiles_per_thread). A cell writes the rows of nodes from x-1 to x+1 (and so in y and z), so
// tiles of one color, at least two cells apart along some axis, also across a
// periodic boundary, never share a row and are added without atomics. Colors
// are taken one after another. All species of a cell are summed into one block
// before it is added, so the rows are written once per cell.

template <ShapeType ST>
void accumulate_species(const Species& species, const Field3d& fieldB, const Domain& domain, const double dt,
                        const std::vector<int>& unknown_index, const ResponseLayout& layout, Operator& response,
                        Field3d* current) {
    static_assert(ST == ShapeType::Linear, "ECSIM response supports Linear shape only");
    RECORD_TIMER;
    std::vector<const ParticlesArray*> charged;
    for (const Species::value_type& entry : species)
        if (!entry.second->is_neutral())
            charged.push_back(entry.second.get());
    if (charged.empty())
        return;

    const Vector3I cells = charged.front()->cells().size();
    for ([[maybe_unused]] const ParticlesArray* sp : charged)
        assert(sp->cells().size() == cells && "species share the cell grid");

    std::array<std::vector<int>, 3> axis_starts;
    Vector3I tiles;
    for (int a = 0; a < 2; ++a) {
        axis_starts[a] = axis_tiles(cells[a], domain.num_cells(a));
        tiles[a] = std::ssize(axis_starts[a]) - 1;
    }

    // z tiles alternate colors too, so a z tile count of 2m multiplies the
    // tiles of one color by m.
    const int xy_tiles_of_color = (tiles.x() + 1) / 2 * ((tiles.y() + 1) / 2);
    const int wanted = min_tiles_per_thread * omp_get_max_threads();
    const int z_pairs = (wanted + xy_tiles_of_color - 1) / xy_tiles_of_color;
    const int interior_z = domain.num_cells(2);
    if (z_pairs <= 1)
        axis_starts[2] = {0, cells.z()};
    else
        axis_starts[2] = axis_tiles(cells.z(), interior_z, std::max(2, interior_z / (2 * z_pairs)));
    tiles[2] = std::ssize(axis_starts[2]) - 1;

    // Tiles by color, three colors per axis.
    std::array<std::vector<Vector3I>, 27> color_tiles;
    for (int tx = 0; tx < tiles.x(); ++tx)
        for (int ty = 0; ty < tiles.y(); ++ty)
            for (int tz = 0; tz < tiles.z(); ++tz) {
                const int color =
                    (tile_color(tx, tiles.x()) * 3 + tile_color(ty, tiles.y())) * 3 + tile_color(tz, tiles.z());
                color_tiles[color].emplace_back(tx, ty, tz);
            }

#pragma omp parallel
    {
        timer::FlatTimer timerOMP("OMP section", -1, timer::MeasureUnit::byte);
        Block tmpBlock;
        PredictedCurrent cellCurrent;
        for (const std::vector<Vector3I>& tiles_of_color : color_tiles) {
#pragma omp for schedule(dynamic, 1)
            for (int t = 0; t < std::ssize(tiles_of_color); ++t) {
                const auto [tx, ty, tz] = tiles_of_color[t].split();
                for (int x = axis_starts[0][tx]; x < axis_starts[0][tx + 1]; ++x)
                    for (int y = axis_starts[1][ty]; y < axis_starts[1][ty + 1]; ++y)
                        for (int z = axis_starts[2][tz]; z < axis_starts[2][tz + 1]; ++z)
                            accumulate_cell(charged, fieldB, domain, dt, unknown_index, layout, Vector3I(x, y, z),
                                            tmpBlock, cellCurrent, response, current);
            }
        }
    }
}

template void accumulate_species<ShapeType::Linear>(const Species&, const Field3d&, const Domain&, double,
                                                    const std::vector<int>&, const ResponseLayout&, Operator&,
                                                    Field3d*);
}   // namespace ecsim_response
