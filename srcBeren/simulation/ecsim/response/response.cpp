#include "response.h"

#include "coupling/stencil.h"
#include "particles/particle_array.h"
#include "response_kernels.h"
#include "runtime/profiler.h"

namespace ecsim_response {
void clear_row_blocks(RowBlocks& row_blocks) {
    // Reuse buffers for an unchanged team size; resize on both growth and shrink.
    if (std::ssize(row_blocks) != omp_get_max_threads())
        row_blocks.resizeAndReset(omp_get_max_threads());
    for (SimpleArrayBuffer<RowBlock<36>>& blocks : row_blocks) blocks.resizeAndReset(0);
}

template <ShapeType ST>
void accumulate_species_reference(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain,
                                  const double dt, BlockMatrix& reference_blocks) {
    RECORD_TIMER;
    static_assert(ST == ShapeType::Linear, "ECSIM response supports Linear shape only");
    if (species.is_neutral())
        return;
#pragma omp parallel for schedule(dynamic, 32)
    for (int pk = 0; pk < species.size(); ++pk)
        for (const Particle& particle : species(pk))
            update_Lmat2_Reference(particle.coord, domain, species.charge(), species.mass(), species.mpw(), fieldB, dt,
                                   reference_blocks);
}

template <typename ColIdx, int DIR, typename Block_t>
static void blockToRowBlocks2Core(int rowIdx, int i_cell, int j_cell, int k_cell, const Block_t& block,
                                  [[maybe_unused]] int Nx, int Ny, int Nz, double tolerance, RowBlock<12>& rowBlock) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };

    for (int x2 = 0; x2 < ColIdx::size_x; ++x2) {
        for (int y2 = 0; y2 < ColIdx::size_y; ++y2) {
            for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                const int colIdx = ColIdx::calculate(x2, y2, z2);
                const double val = block(rowIdx, colIdx, DIR);

                if (std::abs(val) > tolerance) {
                    const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                         k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                    rowBlock.push_back_value(col, val);
                }
            }
        }
    }
}

template <typename RowIdx, int baseDir, typename Block_t>
static void blockToRowBlocks2(int i_cell, int j_cell, int k_cell, const Block_t& block, [[maybe_unused]] int Nx, int Ny,
                              int Nz, double tolerance, SimpleArrayBuffer<RowBlock<36>>& rowBlocks) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };

    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row = vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                     k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                const int rowIdx = RowIdx::calculate(x1, y1, z1);

                RowBlock<12> rowBuffers[3]{row, row, row};
                blockToRowBlocks2Core<XIndexer, baseDir + 0>(rowIdx, i_cell, j_cell, k_cell, block, Nx, Ny, Nz,
                                                             tolerance, rowBuffers[0]);
                blockToRowBlocks2Core<YIndexer, baseDir + 1>(rowIdx, i_cell, j_cell, k_cell, block, Nx, Ny, Nz,
                                                             tolerance, rowBuffers[1]);
                blockToRowBlocks2Core<ZIndexer, baseDir + 2>(rowIdx, i_cell, j_cell, k_cell, block, Nx, Ny, Nz,
                                                             tolerance, rowBuffers[2]);
                const bool isEmpty = rowBuffers[0].nnz == 0 && rowBuffers[1].nnz == 0 && rowBuffers[2].nnz == 0;
                if (!isEmpty) {
                    RowBlock<36> tmp;
                    tmp.mergeFromOthers<3>(rowBuffers);
                    rowBlocks.push_back(tmp);
                }
            }
        }
    }
}

template <ShapeType ST>
void accumulate_species(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, const double dt,
                        RowBlocks& rowBlocksGlobal) {
    static_assert(ST == ShapeType::Linear, "ECSIM response supports Linear shape only");
    if (species.is_neutral())
        return;
    RECORD_TIMER;
#pragma omp parallel
    {
        timer::flatTimer timerOMP("OMP section", -1, timer::MeasureUnit::byte);
        BlockStack tmpBlock;

        SimpleArrayBuffer<RowBlock<36>>& rowBlockThrLocal = rowBlocksGlobal[omp_get_thread_num()];

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

            const Vector3I cell = to_storage_index(cell_origin, domain.grid.ghost_cells());
            const int i = cell.x(), j = cell.y(), k = cell.z();

            constexpr double TOL = 1e-16;

            const int xSize = domain.size().x();
            const int ySize = domain.size().y();
            const int zSize = domain.size().z();

            blockToRowBlocks2<XIndexer, 0>(i, j, k, tmpBlock, xSize, ySize, zSize, TOL, rowBlockThrLocal);
            blockToRowBlocks2<YIndexer, 3>(i, j, k, tmpBlock, xSize, ySize, zSize, TOL, rowBlockThrLocal);
            blockToRowBlocks2<ZIndexer, 6>(i, j, k, tmpBlock, xSize, ySize, zSize, TOL, rowBlockThrLocal);
        }
        timerOMP.m = rowBlockThrLocal.size() * sizeof(rowBlockThrLocal[0]);
    }
}

template void accumulate_species<ShapeType::Linear>(const ParticlesArray&, const Field3d&, const Domain&, double,
                                                    RowBlocks&);
template void accumulate_species_reference<ShapeType::Linear>(const ParticlesArray&, const Field3d&, const Domain&,
                                                              double, BlockMatrix&);
}   // namespace ecsim_response
