#include "response_kernels.h"

#include "coupling/stencil.h"
#include "runtime/profiler.h"

namespace ecsim_response {
void update_Lmat2_Reference(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                            const Field3d& fieldB, const double dt, BlockMatrix& LmatX2) {
    const int SMAX = 2;
    const Vector3R coordLoc = domain.to_cell_coordinates(coord);
    const Stencil3D<2> no = make_point_stencil<ShapeType::Linear>(coordLoc);
    const Stencil3D<2> sh = make_point_stencil<ShapeType::Linear>(coordLoc, Vector3R(0.5));
    // Raw field/block storage is entered only after logical support selection.
    const Vector3I cell = domain.grid.to_storage_index(no.origin);
    const Vector3I cell05 = domain.grid.to_storage_index(sh.origin);
    const int cellLocX = cell.x(), cellLocY = cell.y(), cellLocZ = cell.z();
    const int cellLocX05 = cell05.x(), cellLocY05 = cell05.y(), cellLocZ05 = cell05.z();

    Vector3R B = Vector3R(0.);
    // Keep the reference i/j/k gather order and multiplication order.
    for (int i = 0; i < SMAX; ++i) {
        const int indx = cellLocX + i;
        const int indx05 = cellLocX05 + i;
        for (int j = 0; j < SMAX; ++j) {
            const int indy = cellLocY + j;
            const int indy05 = cellLocY05 + j;
            for (int k = 0; k < SMAX; ++k) {
                const int indz = cellLocZ + k;
                const int indz05 = cellLocZ05 + k;
                const double wx = no(i, Axis::X) * sh(j, Axis::Y) * sh(k, Axis::Z);
                const double wy = sh(i, Axis::X) * no(j, Axis::Y) * sh(k, Axis::Z);
                const double wz = sh(i, Axis::X) * sh(j, Axis::Y) * no(k, Axis::Z);
                B.x() += (wx * fieldB(indx, indy05, indz05, 0));
                B.y() += (wy * fieldB(indx05, indy, indz05, 1));
                B.z() += (wz * fieldB(indx05, indy05, indz, 2));
            }
        }
    }
    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squaredNorm());
    const double betaL = 0.25 * dt * dt * q_m * betaI;

    const int blockIndex = domain.sind(cellLocX, cellLocY, cellLocZ);
    Block& currentBlock = LmatX2[blockIndex];

    const int xOffset = sh.origin.x() - no.origin.x() - XIndexer::offset_x;
    const int yOffset = sh.origin.y() - no.origin.y() - YIndexer::offset_y;
    const int zOffset = sh.origin.z() - no.origin.z() - ZIndexer::offset_z;

    const double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                               {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                               {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    for (int i = 0; i < SMAX; ++i) {
        for (int j = 0; j < SMAX; ++j) {
            for (int k = 0; k < SMAX; ++k) {
                const double s1[3] = {
                    sh(i, Axis::X) * no(j, Axis::Y) * no(k, Axis::Z),   // X
                    no(i, Axis::X) * sh(j, Axis::Y) * no(k, Axis::Z),   // Y
                    no(i, Axis::X) * no(j, Axis::Y) * sh(k, Axis::Z)    // Z
                };
                const int idx1[3] = {BlockDims::indX(xOffset + i, j, k), BlockDims::indY(i, yOffset + j, k),
                                     BlockDims::indZ(i, j, zOffset + k)};

                for (int i1 = 0; i1 < SMAX; ++i1) {
                    for (int j1 = 0; j1 < SMAX; ++j1) {
                        for (int k1 = 0; k1 < SMAX; ++k1) {
                            const double s2[3] = {sh(i1, Axis::X) * no(j1, Axis::Y) * no(k1, Axis::Z),
                                                  no(i1, Axis::X) * sh(j1, Axis::Y) * no(k1, Axis::Z),
                                                  no(i1, Axis::X) * no(j1, Axis::Y) * sh(k1, Axis::Z)};
                            const int idx2[3] = {BlockDims::indX(xOffset + i1, j1, k1),
                                                 BlockDims::indY(i1, yOffset + j1, k1),
                                                 BlockDims::indZ(i1, j1, zOffset + k1)};

                            for (int c1 = 0; c1 < 3; ++c1) {
                                const int rowIndex = idx1[c1];
                                for (int c2 = 0; c2 < 3; ++c2) {
                                    const int colIndex = idx2[c2];
                                    currentBlock(rowIndex, colIndex, c1 * 3 + c2) +=
                                        betaL * s1[c1] * s2[c2] * matB[c1][c2];
                                }
                            }
                        }   // k1
                    }   // j1
                }   // i1
            }   // k
        }   // j
    }   // i
}

void update_Lmat2(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                  const Field3d& fieldB, const double dt, BlockStack& currentBlock) {
    const int SMAX = 2;
    const Vector3R coordLoc = domain.to_cell_coordinates(coord);
    const Stencil3D<2> no = make_point_stencil<ShapeType::Linear>(coordLoc);
    const Stencil3D<2> sh = make_point_stencil<ShapeType::Linear>(coordLoc, Vector3R(0.5));
    // Raw field/block storage is entered only after logical support selection.
    const Vector3I cell = domain.grid.to_storage_index(no.origin);
    const Vector3I cell05 = domain.grid.to_storage_index(sh.origin);
    const int cellLocX = cell.x(), cellLocY = cell.y(), cellLocZ = cell.z();
    const int cellLocX05 = cell05.x(), cellLocY05 = cell05.y(), cellLocZ05 = cell05.z();

    Vector3R B = Vector3R(0.);
    // Keep the reference i/j/k gather order and multiplication order.
    for (int i = 0; i < SMAX; ++i) {
        const int indx = cellLocX + i;
        const int indx05 = cellLocX05 + i;
        for (int j = 0; j < SMAX; ++j) {
            const int indy = cellLocY + j;
            const int indy05 = cellLocY05 + j;
            for (int k = 0; k < SMAX; ++k) {
                const int indz = cellLocZ + k;
                const int indz05 = cellLocZ05 + k;
                const double wx = no(i, Axis::X) * sh(j, Axis::Y) * sh(k, Axis::Z);
                const double wy = sh(i, Axis::X) * no(j, Axis::Y) * sh(k, Axis::Z);
                const double wz = sh(i, Axis::X) * sh(j, Axis::Y) * no(k, Axis::Z);
                B.x() += (wx * fieldB(indx, indy05, indz05, 0));
                B.y() += (wy * fieldB(indx05, indy, indz05, 1));
                B.z() += (wz * fieldB(indx05, indy05, indz, 2));
            }
        }
    }
    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squaredNorm());
    const double betaL = 0.25 * dt * dt * q_m * betaI;

    const int xOffset = sh.origin.x() - no.origin.x() - XIndexer::offset_x;
    const int yOffset = sh.origin.y() - no.origin.y() - YIndexer::offset_y;
    const int zOffset = sh.origin.z() - no.origin.z() - ZIndexer::offset_z;

    double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                         {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                         {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            matB[i][j] *= betaL;
        }
    }

    Vector3d sAll[SMAX * SMAX * SMAX];
    Vector3i idxAll[SMAX * SMAX * SMAX];

    for (int i = 0; i < SMAX; ++i) {
        for (int j = 0; j < SMAX; ++j) {
            for (int k = 0; k < SMAX; ++k) {
                const int ix = (i * SMAX + j) * SMAX + k;
                sAll[ix] = {
                    sh(i, Axis::X) * no(j, Axis::Y) * no(k, Axis::Z),
                    no(i, Axis::X) * sh(j, Axis::Y) * no(k, Axis::Z),
                    no(i, Axis::X) * no(j, Axis::Y) * sh(k, Axis::Z),
                };

                idxAll[ix] = {
                    BlockDims::indX(xOffset + i, j, k),
                    BlockDims::indY(i, yOffset + j, k),
                    BlockDims::indZ(i, j, zOffset + k),
                };
            }
        }
    }
    for (int i1 = 0; i1 < SMAX; ++i1) {
        for (int j1 = 0; j1 < SMAX; ++j1) {
            for (int k1 = 0; k1 < SMAX; ++k1) {
                const int ix1 = (i1 * SMAX + j1) * SMAX + k1;
                const Vector3d& s1 = sAll[ix1];
                const Vector3i& idx1 = idxAll[ix1];

                for (int i2 = 0; i2 < SMAX; ++i2) {
                    for (int j2 = 0; j2 < SMAX; ++j2) {
                        for (int k2 = 0; k2 < SMAX; ++k2) {
                            const int ix2 = (i2 * SMAX + j2) * SMAX + k2;
                            const Vector3d& s2 = sAll[ix2];
                            const Vector3i& idx2 = idxAll[ix2];

                            for (int c1 = 0; c1 < 3; ++c1) {
                                const int rowIndex = idx1[c1];
                                for (int c2 = 0; c2 < 3; ++c2) {
                                    const int colIndex = idx2[c2];
                                    currentBlock(rowIndex, colIndex, c1 * 3 + c2) += s1[c1] * s2[c2] * matB[c1][c2];
                                }
                            }
                        }   // k2
                    }   // j2
                }   // i2
            }   // k1
        }   // j1
    }   // i1
}

/* The same as update_Lmat2, but deals with an array of coordinates => much better CPU utilization
 */
template <int size>
void update_Lmat2(const Vector3R* coord, const Domain& domain, double charge, double mass, double mpw,
                  const Field3d& fieldB, const double dt, BlockStack& currentBlock) {
    constexpr int SMAX = 2;
    static_assert(size > 0);

    const Vector3R firstCoordLoc = domain.to_cell_coordinates(coord[0]);
    const Stencil3D<2> firstNo = make_point_stencil<ShapeType::Linear>(firstCoordLoc);
    const Stencil3D<2> firstSh = make_point_stencil<ShapeType::Linear>(firstCoordLoc, Vector3R(0.5));
    const Vector3I cell = domain.grid.to_storage_index(firstNo.origin);
    const Vector3I cell05 = domain.grid.to_storage_index(firstSh.origin);
    const int cellLocX = cell.x(), cellLocY = cell.y(), cellLocZ = cell.z();
    const int cellLocX05 = cell05.x(), cellLocY05 = cell05.y(), cellLocZ05 = cell05.z();

    Eigen::Matrix<double, size, SMAX> sx;
    Eigen::Matrix<double, size, SMAX> sy;
    Eigen::Matrix<double, size, SMAX> sz;
    Eigen::Matrix<double, size, SMAX> sx05;
    Eigen::Matrix<double, size, SMAX> sy05;
    Eigen::Matrix<double, size, SMAX> sz05;

    for (int i = 0; i < size; ++i) {
        const Vector3R coordLoc = domain.to_cell_coordinates(coord[i]);
        const Stencil3D<2> no = make_point_stencil<ShapeType::Linear>(coordLoc);
        const Stencil3D<2> sh = make_point_stencil<ShapeType::Linear>(coordLoc, Vector3R(0.5));
        // Every batch shares both Yee supports, not only the physical cell.
        assert(no.origin == firstNo.origin);
        assert(sh.origin == firstSh.origin);
        for (int local = 0; local < SMAX; ++local) {
            sx(i, local) = no(local, Axis::X);
            sy(i, local) = no(local, Axis::Y);
            sz(i, local) = no(local, Axis::Z);
            sx05(i, local) = sh(local, Axis::X);
            sy05(i, local) = sh(local, Axis::Y);
            sz05(i, local) = sh(local, Axis::Z);
        }
    }

    std::array<Vector3R, size> B;

    Eigen::Vector<double, SMAX * SMAX * SMAX> xBuf;
    Eigen::Vector<double, SMAX * SMAX * SMAX> yBuf;
    Eigen::Vector<double, SMAX * SMAX * SMAX> zBuf;
    for (int i = 0; i < SMAX; ++i) {
        const int indx = cellLocX + i;
        const int indx05 = cellLocX05 + i;
        for (int j = 0; j < SMAX; ++j) {
            const int indy = cellLocY + j;
            const int indy05 = cellLocY05 + j;
            for (int k = 0; k < SMAX; ++k) {
                const int ix = (i * 2 + j) * 2 + k;
                const int indz = cellLocZ + k;
                const int indz05 = cellLocZ05 + k;
                xBuf[ix] = fieldB(indx, indy05, indz05, 0);
                yBuf[ix] = fieldB(indx05, indy, indz05, 1);
                zBuf[ix] = fieldB(indx05, indy05, indz, 2);
            }
        }
    }

    for (int ix = 0; ix < size; ++ix) {
        B[ix] = Vector3R(0.0);
        for (int i = 0; i < SMAX; ++i) {
            for (int j = 0; j < SMAX; ++j) {
                const double wx0 = sx(ix, i) * sy05(ix, j);
                const double wy0 = sx05(ix, i) * sy(ix, j);
                const double wz0 = sx05(ix, i) * sy05(ix, j);
                for (int k = 0; k < SMAX; ++k) {
                    const int ix2 = (i * 2 + j) * 2 + k;
                    const double tmpX = xBuf[ix2];
                    const double tmpY = yBuf[ix2];
                    const double tmpZ = zBuf[ix2];
                    const double wx = wx0 * sz05(ix, k);
                    const double wy = wy0 * sz05(ix, k);
                    const double wz = wz0 * sz(ix, k);
                    B[ix].x() += wx * tmpX;
                    B[ix].y() += wy * tmpY;
                    B[ix].z() += wz * tmpZ;
                }
            }
        }
    }

    const double q_m = charge / mass;
    const int xOffset = firstSh.origin.x() - firstNo.origin.x() - XIndexer::offset_x;
    const int yOffset = firstSh.origin.y() - firstNo.origin.y() - YIndexer::offset_y;
    const int zOffset = firstSh.origin.z() - firstNo.origin.z() - ZIndexer::offset_z;

    Eigen::Matrix<Eigen::Vector<double, size>, 3, 3> matB;
    for (int i = 0; i < size; ++i) {
        const Vector3R b = 0.5 * dt * q_m * B[i];
        const double betaI = mpw * charge / (1.0 + b.squaredNorm());
        const double betaL = 0.25 * dt * dt * q_m * betaI;
        double tmp[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                            {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                            {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                matB(j, k)(i) = betaL * tmp[j][k];
            }
        }
    }

    Eigen::Vector<Eigen::Matrix<double, size, 3>, SMAX * SMAX * SMAX> sAll;

    Vector3i idxAll[SMAX * SMAX * SMAX];

    for (int i = 0; i < SMAX; ++i) {
        for (int j = 0; j < SMAX; ++j) {
            for (int k = 0; k < SMAX; ++k) {
                const int ix = (i * SMAX + j) * SMAX + k;
                for (int l = 0; l < size; ++l) {
                    sAll(ix).row(l) = Eigen::Vector3d{
                        sx05(l, i) * sy(l, j) * sz(l, k),
                        sx(l, i) * sy05(l, j) * sz(l, k),
                        sx(l, i) * sy(l, j) * sz05(l, k),
                    };
                }

                idxAll[ix] = {
                    BlockDims::indX(xOffset + i, j, k),
                    BlockDims::indY(i, yOffset + j, k),
                    BlockDims::indZ(i, j, zOffset + k),
                };
            }
        }
    }

    for (int i1 = 0; i1 < SMAX; ++i1) {
        for (int j1 = 0; j1 < SMAX; ++j1) {
            for (int k1 = 0; k1 < SMAX; ++k1) {
                const int ix1 = (i1 * SMAX + j1) * SMAX + k1;
                const Vector3i& idx1 = idxAll[ix1];
                const Eigen::Matrix<double, size, 3>& sLoc1 = sAll(ix1);

                for (int i2 = 0; i2 < SMAX; ++i2) {
                    for (int j2 = 0; j2 < SMAX; ++j2) {
                        for (int k2 = 0; k2 < SMAX; ++k2) {
                            const int ix2 = (i2 * SMAX + j2) * SMAX + k2;
                            const Vector3i& idx2 = idxAll[ix2];
                            const Eigen::Matrix<double, size, 3>& sLoc2 = sAll(ix2);

                            for (int c1 = 0; c1 < 3; ++c1) {
                                const int rowIndex = idx1[c1];
                                for (int c2 = 0; c2 < 3; ++c2) {
                                    const int colIndex = idx2[c2];
                                    const Eigen::Vector<double, size>& matBLoc = matB(c1, c2);
                                    double acc = 0.0;
#pragma omp simd reduction(+ : acc)
                                    for (int ix = 0; ix < size; ix += 1) {
                                        acc += sLoc1(ix, c1) * sLoc2(ix, c2) * matBLoc(ix);
                                    }
                                    currentBlock(rowIndex, colIndex, c1 * 3 + c2) += acc;
                                }
                            }
                        }   // k2
                    }   // j2
                }   // i2
            }   // k1
        }   // j1
    }   // i1
}

template void update_Lmat2<2>(const Vector3R* coord, const Domain& domain, double charge, double mass, double mpw,
                              const Field3d& fieldB, const double dt, BlockStack& currentBlock);
template void update_Lmat2<4>(const Vector3R* coord, const Domain& domain, double charge, double mass, double mpw,
                              const Field3d& fieldB, const double dt, BlockStack& currentBlock);
template void update_Lmat2<8>(const Vector3R* coord, const Domain& domain, double charge, double mass, double mpw,
                              const Field3d& fieldB, const double dt, BlockStack& currentBlock);
template void update_Lmat2<16>(const Vector3R* coord, const Domain& domain, double charge, double mass, double mpw,
                               const Field3d& fieldB, const double dt, BlockStack& currentBlock);

}   // namespace ecsim_response
