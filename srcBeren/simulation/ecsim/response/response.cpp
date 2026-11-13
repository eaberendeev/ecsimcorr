#include "response.h"

#include <algorithm>

#include "coupling/stencil.h"
#include "particles/particle_array.h"
#include "response_kernels.h"
#include "runtime/profiler.h"

namespace ecsim_response {
Operator response_pattern(const Domain& domain, const std::array<bool, 3>& periodic,
                          const std::vector<int>& unknown_index, const std::vector<int>& storage_index) {
    RECORD_TIMER;
    const Grid& grid = domain.grid;
    const int n = std::ssize(storage_index);
    std::vector<std::vector<int>> columns(n);
#pragma omp parallel for schedule(dynamic, 1024)
    for (int row = 0; row < n; ++row) {
        const std::array<int, 4> node = grid.pos_vind_range(storage_index[row]);
        const int c1 = node[3];
        for (int c2 = 0; c2 < 3; ++c2) {
            Vector3I lo, hi;
            for (int a = 0; a < 3; ++a) {
                lo[a] = (a == c2 && a != c1) ? -2 : -1;
                hi[a] = (a == c1 && a != c2) ? 2 : 1;
            }
            for (int dx = lo.x(); dx <= hi.x(); ++dx)
                for (int dy = lo.y(); dy <= hi.y(); ++dy)
                    for (int dz = lo.z(); dz <= hi.z(); ++dz) {
                        Vector3I neighbor(node[0] + dx, node[1] + dy, node[2] + dz);
                        bool stored = true;
                        for (int a = 0; a < 3; ++a) {
                            if (periodic[a])
                                neighbor[a] = grid.periodic_node_owner(neighbor[a], a);
                            stored = stored && neighbor[a] >= 0 && neighbor[a] < grid.size(a);
                        }
                        const int column =
                            stored ? unknown_index[domain.vind(neighbor.x(), neighbor.y(), neighbor.z(), c2)] : -1;
                        if (column >= 0)
                            columns[row].push_back(column);
                    }
        }
        std::sort(columns[row].begin(), columns[row].end());
        columns[row].erase(std::unique(columns[row].begin(), columns[row].end()), columns[row].end());
    }
    std::vector<Trip> trips;
    for (int row = 0; row < n; ++row)
        for (int column : columns[row]) trips.emplace_back(row, column, 0.0);
    Operator pattern(n, n);
    pattern.setFromTriplets(trips.begin(), trips.end());
    pattern.makeCompressed();
    return pattern;
}

// Add a sorted row of contributions to the same row of `response`.
static void add_row(const RowBlock<36>& block, Operator& response) {
    const int begin = response.outerIndexPtr()[block.row];
    const int* columns = response.innerIndexPtr() + begin;
    double* values = response.valuePtr() + begin;
    int p = 0;
    for (int e = 0; e < block.nnz; ++e) {
        while (columns[p] < block.columns[e]) ++p;
        assert(columns[p] == block.columns[e] && "response entry outside the fixed pattern");
#pragma omp atomic update
        values[p] += block.values[e];
    }
}

// Columns are field unknowns: periodic images of a node map to the unknown of its owner,
// which can break the column order and join columns on short periodic axes. Fixed
// components (outside the field domain, index -1) carry no response.
template <typename ColIdx, int DIR, typename Block_t>
static void collect_block_row(int rowIdx, int i_cell, int j_cell, int k_cell, const Block_t& block,
                              const Domain& domain, const std::vector<int>& unknown_index, double tolerance,
                              RowBlock<12>& rowBlock) {
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
        rowBlock.sort_and_combine();
}

template <typename RowIdx, int baseDir, typename Block_t>
static void add_block_rows(int i_cell, int j_cell, int k_cell, const Block_t& block, const Domain& domain,
                           const std::vector<int>& unknown_index, double tolerance, Operator& response) {
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row =
                    unknown_index[domain.vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                              k_cell + z1 + RowIdx::offset_z, RowIdx::dir)];

                if (row < 0)
                    continue;
                const int rowIdx = RowIdx::calculate(x1, y1, z1);

                RowBlock<12> rowBuffers[3]{row, row, row};
                collect_block_row<XIndexer, baseDir + 0>(rowIdx, i_cell, j_cell, k_cell, block, domain, unknown_index,
                                                         tolerance, rowBuffers[0]);
                collect_block_row<YIndexer, baseDir + 1>(rowIdx, i_cell, j_cell, k_cell, block, domain, unknown_index,
                                                         tolerance, rowBuffers[1]);
                collect_block_row<ZIndexer, baseDir + 2>(rowIdx, i_cell, j_cell, k_cell, block, domain, unknown_index,
                                                         tolerance, rowBuffers[2]);
                const bool isEmpty = rowBuffers[0].nnz == 0 && rowBuffers[1].nnz == 0 && rowBuffers[2].nnz == 0;
                if (!isEmpty) {
                    RowBlock<36> tmp;
                    tmp.mergeFromOthers<3>(rowBuffers);
                    add_row(tmp, response);
                }
            }
        }
    }
}

template <ShapeType ST>
void accumulate_species(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, const double dt,
                        const std::vector<int>& unknown_index, Operator& response) {
    static_assert(ST == ShapeType::Linear, "ECSIM response supports Linear shape only");
    if (species.is_neutral())
        return;
    RECORD_TIMER;
#pragma omp parallel
    {
        Block tmpBlock;

#pragma omp for schedule(dynamic, species.omp_granularity())
        for (auto pk = 0; pk < species.size(); ++pk) {
            const std::vector<Particle>& currVec = species(pk);
            if (currVec.size() == 0) {
                continue;
            }
            tmpBlock.setZero();
            const Vector3R first_coord = domain.to_cell_coordinates(currVec.front().coord);
            const Vector3I cell_origin = checked_floor_index(first_coord);

            if (currVec.size() >= 64) {
                constexpr int buffSize = 16;
                std::array<std::array<Vector3R, buffSize>, 8> coordBuffers;
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

                    coordBuffers[currKind][accumulatedCount[currKind]] = coord;
                    accumulatedCount[currKind] += 1;

                    if (accumulatedCount[currKind] == buffSize) {
                        update_Lmat2<buffSize>(&coordBuffers[currKind][0], domain, species.charge(), species.mass(),
                                               species.mpw(), fieldB, dt, tmpBlock);
                        accumulatedCount[currKind] = 0;
                    }
                }

                for (int kind = 0; kind < 8; ++kind) {
                    const int accumulated = accumulatedCount[kind];
                    int offset = 0;
                    if (accumulated >= 8) {
                        update_Lmat2<8>(&coordBuffers[kind][offset], domain, species.charge(), species.mass(),
                                        species.mpw(), fieldB, dt, tmpBlock);
                        offset += 8;
                    }
                    if (accumulated - offset >= 4) {
                        update_Lmat2<4>(&coordBuffers[kind][offset], domain, species.charge(), species.mass(),
                                        species.mpw(), fieldB, dt, tmpBlock);
                        offset += 4;
                    }
                    if (accumulated - offset >= 2) {
                        update_Lmat2<2>(&coordBuffers[kind][offset], domain, species.charge(), species.mass(),
                                        species.mpw(), fieldB, dt, tmpBlock);
                        offset += 2;
                    }
                    if (accumulated - offset >= 1) {
                        update_Lmat2(coordBuffers[kind][offset], domain, species.charge(), species.mass(),
                                     species.mpw(), fieldB, dt, tmpBlock);
                        offset += 1;
                    }
                    assert(accumulated == offset);
                }
            } else {
                for (const Particle& particle : currVec) {
                    const Vector3R& coord = particle.coord;
                    update_Lmat2(coord, domain, species.charge(), species.mass(), species.mpw(), fieldB, dt, tmpBlock);
                }
            }

            const Vector3I cell = to_storage_index(cell_origin);
            const int i = cell.x(), j = cell.y(), k = cell.z();

            constexpr double TOL = 1e-16;

            add_block_rows<XIndexer, 0>(i, j, k, tmpBlock, domain, unknown_index, TOL, response);
            add_block_rows<YIndexer, 3>(i, j, k, tmpBlock, domain, unknown_index, TOL, response);
            add_block_rows<ZIndexer, 6>(i, j, k, tmpBlock, domain, unknown_index, TOL, response);
        }
    }
}

template void accumulate_species<ShapeType::Linear>(const ParticlesArray&, const Field3d&, const Domain&, double,
                                                    const std::vector<int>&, Operator&);
}   // namespace ecsim_response
