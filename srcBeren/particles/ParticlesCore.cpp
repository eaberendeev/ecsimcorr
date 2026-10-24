#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "containers.h"
#include "memory.h"

void ParticlesArray::move(double dt) {
    timer::commonTimer timer(std::source_location::current().function_name(), -1, timer::MeasureUnit::byte);
    int64_t totalParticles = 0;
#pragma omp parallel for schedule(dynamic, 32) reduction(+ : totalParticles)
    for (auto k = 0; k < size(); ++k) {
        totalParticles += std::ssize(particlesData(k));
        for (auto& particle : particlesData(k)) {
            particle.move(dt);
        }
    }
    timer.flat.m = sizeof(Particle) * totalParticles;
}

// Fill Lmatrix by each particles
template <ShapeType ST>
void ParticlesArray::fill_matrixL2_Reference(Mesh& mesh, const Field3d& fieldB, const Domain& domain,
                                             const double dt) const {
    RECORD_TIMER;
    if (is_neutral())
        return;

    static_assert(ST == ShapeType::Linear, "fill_matrixL2_Reference supports Linear shape only");
    fill_matrixL_impl_linear2_Reference(mesh, fieldB, domain, dt);
}

// Fill Lmatrix by each particles, produce row-block for each cell
template <ShapeType ST>
void ParticlesArray::fill_matrixL2(Mesh& mesh, const Field3d& fieldB, const Domain& domain, const double dt,
                                   SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>& rowBlocksGlobal) const {
    RECORD_TIMER;
    if (is_neutral())
        return;

    static_assert(ST == ShapeType::Linear, "fill_matrixL2 is designed for now only for Linear shape");
    fill_matrixL_impl_linear2(mesh, fieldB, domain, dt, rowBlocksGlobal);
}

template void ParticlesArray::fill_matrixL2_Reference<ShapeType::Linear>(Mesh&, const Field3d&, const Domain&,
                                                                        const double) const;
template void ParticlesArray::fill_matrixL2<ShapeType::Linear>(Mesh&, const Field3d&, const Domain&, const double,
                                                              SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>&) const;

void ParticlesArray::fill_matrixL_impl_linear2_Reference(Mesh& mesh, const Field3d& fieldB, const Domain& domain,
                                                         const double dt) const {
    RECORD_TIMER;

#pragma omp parallel for schedule(dynamic, 32)
    for (auto pk = 0; pk < size(); ++pk) {
        for (auto& particle : particlesData(pk)) {
            const auto coord = particle.coord;
            mesh.update_Lmat2_Reference(coord, domain, charge, mass_, mpw_, fieldB, dt);
        }
    }
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

void ParticlesArray::fill_matrixL_impl_linear2(
    const Mesh& mesh, const Field3d& fieldB, const Domain& domain, const double dt,
    SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>& rowBlocksGlobal) const {
    RECORD_TIMER;
#pragma omp parallel
    {
        timer::flatTimer timerOMP("OMP section", -1, timer::MeasureUnit::byte);
        BlockStack tmpBlock;

        SimpleArrayBuffer<RowBlock<36>>& rowBlockThrLocal = rowBlocksGlobal[omp_get_thread_num()];

        // Distribute cells in small chunks to keep threads busy on narrow domains.
#pragma omp for schedule(dynamic, 8)
        for (auto pk = 0; pk < size(); ++pk) {
            const std::vector<Particle>& currVec = particlesData(pk);
            if (currVec.size() == 0) {
                continue;
            }

            tmpBlock.setZero();
            for (const auto& particle : currVec) {
                mesh.update_Lmat2(particle.coord, domain, charge, mass_, mpw_, fieldB, dt, tmpBlock);
            }

            const Vector3R coord = currVec[0].coord;
            const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
            const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
            const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;

            const int i = int(coordLocX);
            const int j = int(coordLocY);
            const int k = int(coordLocZ);

            constexpr double TOL = 1e-16;

            const int xSize = mesh.sizes().x();
            const int ySize = mesh.sizes().y();
            const int zSize = mesh.sizes().z();

            // X component
            blockToRowBlocks2<XIndexer, 0>(i, j, k, tmpBlock, xSize, ySize, zSize, TOL, rowBlockThrLocal);
            blockToRowBlocks2<YIndexer, 3>(i, j, k, tmpBlock, xSize, ySize, zSize, TOL, rowBlockThrLocal);
            blockToRowBlocks2<ZIndexer, 6>(i, j, k, tmpBlock, xSize, ySize, zSize, TOL, rowBlockThrLocal);
        }
        timerOMP.m = rowBlockThrLocal.size() * sizeof(rowBlockThrLocal[0]);
    }
}
