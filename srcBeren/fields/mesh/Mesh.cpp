#include "Mesh.h"

#include "Shape.h"
#include "World.h"
#include "interpolation.h"
#include "row_block.h"
#include "solverSLE.h"
#include "thread_partitioned_matrix.h"
#include "timer.h"

void Mesh::init(const Domain& domain, double dt, BoundaryConditionHandler& bc_handler) {
    RECORD_TIMER;

    Lmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Lmat2.resize(domain.total_size() * 3, domain.total_size() * 3);
    Mmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Imat.resize(domain.total_size() * 3, domain.total_size() * 3);
    IMmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    chargeDensityOld.resize(domain.size(), 1);
    chargeDensity.resize(domain.size(), 1);
    divE.resize(domain.total_size(), domain.total_size() * 3);
    // TODO: move sind, converter func to BlockMatrix, resize with 3dim
    LmatX2.resize(domain.total_size());
    LmatX_NGP.resize(domain.total_size());

    xCellSize = domain.cell_size().x();
    yCellSize = domain.cell_size().y();
    zCellSize = domain.cell_size().z();
    xSize = domain.size().x();
    ySize = domain.size().y();
    zSize = domain.size().z();

    Operator curlBtmp;
    Operator curlEtmp;
    curlBtmp.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlEtmp.resize(domain.total_size() * 3, domain.total_size() * 3);

    stencil_Imat(Imat, domain);
    stencil_curlE(curlEtmp, domain, bc_handler);
    stencil_curlB(curlBtmp, domain, bc_handler);

    curlE = ThreadPartitionedSparseMatrix<double>(curlEtmp);
    curlB = ThreadPartitionedSparseMatrix<double>(curlBtmp);

    stencil_divE(divE, domain, bc_handler);

    Mmat = -0.25 * dt * dt * curlBtmp * curlEtmp;
    IMmat = Imat - Mmat;
    IMmat.makeCompressed();
}

void Mesh::print_operator(const Operator& oper) {
    for (int k = 0; k < oper.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(oper, k); it; ++it) {
            std::cout << pos_vind(it.row(), 0) << " " << pos_vind(it.row(), 1) << " " << pos_vind(it.row(), 2) << " "
                      << pos_vind(it.row(), 3) << " " << pos_vind(it.col(), 0) << " " << pos_vind(it.col(), 1) << " "
                      << pos_vind(it.col(), 2) << " " << pos_vind(it.col(), 3) << " " << it.value() << "\n";
        }
    }
}

void Mesh::prepare() {
}

void Mesh::fdtd_explicit(Field3d& E, Field3d& B, const Field3d& J, const double dt) {
    RECORD_TIMER;
    E += 0.5 * dt * (curlB * B) - 0.5 * dt * J;
    B -= 0.5 * dt * (curlE * E);
}

void Mesh::computeB(const Field3d& fieldE, const Field3d& fieldEn, Field3d& fieldB, double dt) {
    RECORD_TIMER;
    fieldB -= (0.5 * dt) * (curlE * (fieldE + fieldEn));
}

void Mesh::compute_fieldB(Field3d& Bn, const Field3d& B, const Field3d& E, const Field3d& En, double dt) {
    RECORD_TIMER;
    Bn = B - (0.5 * dt) * (curlE * (E + En));
}

void Mesh::update_Lmat2_Reference(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                                  const Field3d& fieldB, const double dt) {
    const int SMAX = 2;   // SHAPE_SIZE;
    alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
    alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];

    const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    const double coordLocX05 = coordLocX - 0.5;
    const double coordLocY05 = coordLocY - 0.5;
    const double coordLocZ05 = coordLocZ - 0.5;

    const int cellLocX = int(coordLocX);
    const int cellLocY = int(coordLocY);
    const int cellLocZ = int(coordLocZ);
    const int cellLocX05 = int(coordLocX05);
    const int cellLocY05 = int(coordLocY05);
    const int cellLocZ05 = int(coordLocZ05);

    sx[1] = (coordLocX - cellLocX);
    sx[0] = 1 - sx[1];
    sy[1] = (coordLocY - cellLocY);
    sy[0] = 1 - sy[1];
    sz[1] = (coordLocZ - cellLocZ);
    sz[0] = 1 - sz[1];

    sx05[1] = (coordLocX05 - cellLocX05);
    sx05[0] = 1 - sx05[1];
    sy05[1] = (coordLocY05 - cellLocY05);
    sy05[0] = 1 - sy05[1];
    sz05[1] = (coordLocZ05 - cellLocZ05);
    sz05[0] = 1 - sz05[1];

    Vector3R B = Vector3R(0.);
    // TODO: change to interpolation function
    for (int i = 0; i < SMAX; ++i) {
        const int indx = cellLocX + i;
        const int indx05 = cellLocX05 + i;
        for (int j = 0; j < SMAX; ++j) {
            const int indy = cellLocY + j;
            const int indy05 = cellLocY05 + j;
            for (int k = 0; k < SMAX; ++k) {
                const int indz = cellLocZ + k;
                const int indz05 = cellLocZ05 + k;
                const double wx = sx[i] * sy05[j] * sz05[k];
                const double wy = sx05[i] * sy[j] * sz05[k];
                const double wz = sx05[i] * sy05[j] * sz[k];
                B.x() += (wx * fieldB(indx, indy05, indz05, 0));
                B.y() += (wy * fieldB(indx05, indy, indz05, 1));
                B.z() += (wz * fieldB(indx05, indy05, indz, 2));
            }
        }
    }
    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.25 * dt * dt * q_m * betaI;

    const int blockIndex = sind(cellLocX, cellLocY, cellLocZ);
    Block& currentBlock = LmatX2[blockIndex];

    const int xOffset = cellLocX05 - cellLocX + 1;
    const int yOffset = cellLocY05 - cellLocY + 1;
    const int zOffset = cellLocZ05 - cellLocZ + 1;

    const double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                               {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                               {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    for (int i = 0; i < SMAX; ++i) {
        for (int j = 0; j < SMAX; ++j) {
            for (int k = 0; k < SMAX; ++k) {
                const double s1[3] = {
                    sx05[i] * sy[j] * sz[k],   // X
                    sx[i] * sy05[j] * sz[k],   // Y
                    sx[i] * sy[j] * sz05[k]    // Z
                };
                const int idx1[3] = {BlockDims::indX(xOffset + i, j, k), BlockDims::indY(i, yOffset + j, k),
                                     BlockDims::indZ(i, j, zOffset + k)};

                for (int i1 = 0; i1 < SMAX; ++i1) {
                    for (int j1 = 0; j1 < SMAX; ++j1) {
                        for (int k1 = 0; k1 < SMAX; ++k1) {
                            const double s2[3] = {sx05[i1] * sy[j1] * sz[k1], sx[i1] * sy05[j1] * sz[k1],
                                                  sx[i1] * sy[j1] * sz05[k1]};
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

void Mesh::update_Lmat2(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                        const Field3d& fieldB, const double dt, BlockStack& currentBlock) const {
    const int SMAX = 2;   // SHAPE_SIZE;
    alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
    alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];

    const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    const double coordLocX05 = coordLocX - 0.5;
    const double coordLocY05 = coordLocY - 0.5;
    const double coordLocZ05 = coordLocZ - 0.5;

    const int cellLocX = int(coordLocX);
    const int cellLocY = int(coordLocY);
    const int cellLocZ = int(coordLocZ);
    const int cellLocX05 = int(coordLocX05);
    const int cellLocY05 = int(coordLocY05);
    const int cellLocZ05 = int(coordLocZ05);

    sx[1] = (coordLocX - cellLocX);
    sx[0] = 1 - sx[1];
    sy[1] = (coordLocY - cellLocY);
    sy[0] = 1 - sy[1];
    sz[1] = (coordLocZ - cellLocZ);
    sz[0] = 1 - sz[1];

    sx05[1] = (coordLocX05 - cellLocX05);
    sx05[0] = 1 - sx05[1];
    sy05[1] = (coordLocY05 - cellLocY05);
    sy05[0] = 1 - sy05[1];
    sz05[1] = (coordLocZ05 - cellLocZ05);
    sz05[0] = 1 - sz05[1];

    // timer::flatTimer timerPrelim("preliminary");

    Vector3R B = Vector3R(0.);
    // TODO: change to interpolation function
    for (int i = 0; i < SMAX; ++i) {
        const int indx = cellLocX + i;
        const int indx05 = cellLocX05 + i;
        for (int j = 0; j < SMAX; ++j) {
            const int indy = cellLocY + j;
            const int indy05 = cellLocY05 + j;
            for (int k = 0; k < SMAX; ++k) {
                const int indz = cellLocZ + k;
                const int indz05 = cellLocZ05 + k;
                const double wx = sx[i] * sy05[j] * sz05[k];
                const double wy = sx05[i] * sy[j] * sz05[k];
                const double wz = sx05[i] * sy05[j] * sz[k];
                B.x() += (wx * fieldB(indx, indy05, indz05, 0));
                B.y() += (wy * fieldB(indx05, indy, indz05, 1));
                B.z() += (wz * fieldB(indx05, indy05, indz, 2));
            }
        }
    }
    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.25 * dt * dt * q_m * betaI;

    const int xOffset = cellLocX05 - cellLocX + 1;
    const int yOffset = cellLocY05 - cellLocY + 1;
    const int zOffset = cellLocZ05 - cellLocZ + 1;

    double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                         {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                         {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    Eigen::Matrix3d matBEig;
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
                    sx05[i] * sy[j] * sz[k],
                    sx[i] * sy05[j] * sz[k],
                    sx[i] * sy[j] * sz05[k],
                };

                idxAll[ix] = {
                    BlockDims::indX(xOffset + i, j, k),
                    BlockDims::indY(i, yOffset + j, k),
                    BlockDims::indZ(i, j, zOffset + k),
                };
            }
        }
    }
    // thread_local static int timerCounter = 0;
    // timer::flatTimer timerRest(timer::NoStart{});
    // if (timerCounter < 10000) {
    //     timerRest.start("rest loop ref");
    //     timerCounter += 1;
    // }

    // timerPrelim.finish();
    // timer::flatTimer timerRest("rest loop");
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

template <int maxSize>
void Mesh::update_Lmat2(const Vector3R* coord, int, const Domain& domain, double charge, double mass, double mpw,
                        const Field3d& fieldB, const double dt, BlockStack& currentBlock) const {
    static constexpr int size = maxSize;
    constexpr int SMAX = 2;   // SHAPE_SIZE;
    assert(size >= 0 && size <= maxSize);

    const double firstCoordLocX = coord[0].x() / domain.cell_size().x() + GHOST_CELLS;
    const double firstCoordLocY = coord[0].y() / domain.cell_size().y() + GHOST_CELLS;
    const double firstCoordLocZ = coord[0].z() / domain.cell_size().z() + GHOST_CELLS;
    const double firstCoordLocX05 = firstCoordLocX - 0.5;
    const double firstCoordLocY05 = firstCoordLocY - 0.5;
    const double firstCoordLocZ05 = firstCoordLocZ - 0.5;

    const int cellLocX = int(firstCoordLocX);
    const int cellLocY = int(firstCoordLocY);
    const int cellLocZ = int(firstCoordLocZ);
    const int cellLocX05 = int(firstCoordLocX05);
    const int cellLocY05 = int(firstCoordLocY05);
    const int cellLocZ05 = int(firstCoordLocZ05);

    const int cellLocXd = static_cast<double>(cellLocX);
    const int cellLocYd = static_cast<double>(cellLocY);
    const int cellLocZd = static_cast<double>(cellLocZ);
    const int cellLocXd05 = static_cast<double>(cellLocX05);
    const int cellLocYd05 = static_cast<double>(cellLocY05);
    const int cellLocZd05 = static_cast<double>(cellLocZ05);

    // std::array<Eigen::Vector<double, SMAX>, maxSize> sx;
    // std::array<Eigen::Vector<double, SMAX>, maxSize> sy;
    // std::array<Eigen::Vector<double, SMAX>, maxSize> sz;
    // std::array<Eigen::Vector<double, SMAX>, maxSize> sx05;
    // std::array<Eigen::Vector<double, SMAX>, maxSize> sy05;
    // std::array<Eigen::Vector<double, SMAX>, maxSize> sz05;

    Eigen::Matrix<double, maxSize, SMAX> sx;
    Eigen::Matrix<double, maxSize, SMAX> sy;
    Eigen::Matrix<double, maxSize, SMAX> sz;
    Eigen::Matrix<double, maxSize, SMAX> sx05;
    Eigen::Matrix<double, maxSize, SMAX> sy05;
    Eigen::Matrix<double, maxSize, SMAX> sz05;

    for (int i = 0; i < size; ++i) {
        const double coordLocX = coord[i].x() / domain.cell_size().x() + GHOST_CELLS;
        const double coordLocY = coord[i].y() / domain.cell_size().y() + GHOST_CELLS;
        const double coordLocZ = coord[i].z() / domain.cell_size().z() + GHOST_CELLS;
        const double coordLocX05 = coordLocX - 0.5;
        const double coordLocY05 = coordLocY - 0.5;
        const double coordLocZ05 = coordLocZ - 0.5;

        assert(int(coordLocX) == cellLocX);
        assert(int(coordLocY) == cellLocY);
        assert(int(coordLocZ) == cellLocZ);

        assert(int(coordLocX05) == cellLocX05);
        assert(int(coordLocY05) == cellLocY05);
        assert(int(coordLocZ05) == cellLocZ05);

        sx(i, 1) = (coordLocX - cellLocXd);
        sx(i, 0) = 1 - sx(i, 1);
        sy(i, 1) = (coordLocY - cellLocYd);
        sy(i, 0) = 1 - sy(i, 1);
        sz(i, 1) = (coordLocZ - cellLocZd);
        sz(i, 0) = 1 - sz(i, 1);

        sx05(i, 1) = (coordLocX05 - cellLocXd05);
        sx05(i, 0) = 1 - sx05(i, 1);
        sy05(i, 1) = (coordLocY05 - cellLocYd05);
        sy05(i, 0) = 1 - sy05(i, 1);
        sz05(i, 1) = (coordLocZ05 - cellLocZd05);
        sz05(i, 0) = 1 - sz05(i, 1);
    }

    // timer::flatTimer timerPrelim("preliminary");
    std::array<Vector3R, maxSize> B;
    for (int i = 0; i < size; ++i) {
        B[i] = Vector3R(0.0);
    }

    for (int i = 0; i < SMAX; ++i) {
        const int indx = cellLocX + i;
        const int indx05 = cellLocX05 + i;
        for (int j = 0; j < SMAX; ++j) {
            const int indy = cellLocY + j;
            const int indy05 = cellLocY05 + j;
            for (int k = 0; k < SMAX; ++k) {
                const int indz = cellLocZ + k;
                const int indz05 = cellLocZ05 + k;
                const double tmpX = fieldB(indx, indy05, indz05, 0);
                const double tmpY = fieldB(indx05, indy, indz05, 1);
                const double tmpZ = fieldB(indx05, indy05, indz, 2);
                for (int ix = 0; ix < maxSize; ++ix) {
                    const double wx = sx(ix, i) * sy05(ix, j) * sz05(ix, k);
                    const double wy = sx05(ix, i) * sy(ix, j) * sz05(ix, k);
                    const double wz = sx05(ix, i) * sy05(ix, j) * sz(ix, k);
                    B[ix].x() += wx * tmpX;
                    B[ix].y() += wy * tmpY;
                    B[ix].z() += wz * tmpZ;
                }
            }
        }
    }

    const double q_m = charge / mass;
    const int xOffset = cellLocX05 - cellLocX + 1;
    const int yOffset = cellLocY05 - cellLocY + 1;
    const int zOffset = cellLocZ05 - cellLocZ + 1;

    Eigen::Matrix<Eigen::Vector<double, maxSize>, 3, 3> matB;
    // std::array<double[3][3], maxSize> matB;
    for (int i = 0; i < size; ++i) {
        const Vector3R b = 0.5 * dt * q_m * B[i];
        const double betaI = mpw * charge / (1.0 + b.squared());
        const double betaL = 0.25 * dt * dt * q_m * betaI;
        double tmp[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                            {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                            {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                matB(j, k)(i) = betaL * tmp[j][k];
                // matB[i][j][k] = betaL * tmp[j][k];
            }
        }
    }

    // std::array<Vector3d[SMAX * SMAX * SMAX], maxSize> sAll;

    Eigen::Vector<Eigen::Matrix<double, maxSize, 3>, SMAX * SMAX * SMAX> sAll;
    // std::mdspan<double, std::extents<int, 3, maxSize, SMAX*SMAX*SMAX>, std::layout_left> sAllMdspan(sAllBuf)
    // Eigen::MatrixXd<double, maxSize, SMAX*SMAX*SMAX> sAll;

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
    // timerPrelim.finish();
    // thread_local static int timerCounter = 0;
    // timer::flatTimer timerRest(timer::NoStart{});
    // if (timerCounter < 10000) {
    //     timerRest.start("rest loop test");
    //     timerCounter += 1;
    // }
    for (int i1 = 0; i1 < SMAX; ++i1) {
        for (int j1 = 0; j1 < SMAX; ++j1) {
            for (int k1 = 0; k1 < SMAX; ++k1) {
                const int ix1 = (i1 * SMAX + j1) * SMAX + k1;
                // const Vector3d& s1 = sAll[ix1];
                const Vector3i& idx1 = idxAll[ix1];
                const Eigen::Matrix<double, maxSize, 3>& sLoc1 = sAll(ix1);

                for (int i2 = 0; i2 < SMAX; ++i2) {
                    for (int j2 = 0; j2 < SMAX; ++j2) {
                        for (int k2 = 0; k2 < SMAX; ++k2) {
                            const int ix2 = (i2 * SMAX + j2) * SMAX + k2;
                            // const Vector3d& s2 = sAll[ix2];
                            const Vector3i& idx2 = idxAll[ix2];

                            const Eigen::Matrix<double, maxSize, 3>& sLoc2 = sAll(ix2);

                            for (int c1 = 0; c1 < 3; ++c1) {
                                const int rowIndex = idx1[c1];
                                for (int c2 = 0; c2 < 3; ++c2) {
                                    const int colIndex = idx2[c2];
                                    const Eigen::Vector<double, maxSize>& matBLoc = matB(c1, c2);
                                    double acc1 = 0.0;
                                    double acc2 = 0.0;
                                    static_assert(size % 2 == 0);
#pragma omp simd
                                    for (int ix = 0; ix < size / 2; ix += 1) {
                                        acc1 += sLoc1(ix, c1) * sLoc2(ix, c2) * matBLoc(ix);
                                        acc2 += sLoc1(ix + size / 2, c1) * sLoc2(ix + size / 2, c2) *
                                                matBLoc(ix + size / 2);
                                    }
                                    currentBlock(rowIndex, colIndex, c1 * 3 + c2) += acc1 + acc2;
                                }
                            }
                        }   // k2
                    }   // j2
                }   // i2
            }   // k1
        }   // j1
    }   // i1
}

void Mesh::update_Lmat2_NGP(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                            const Field3d& fieldB, const double dt) {
    RECORD_TIMER;

    Vector3R B = Vector3R(0.);
    const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    const double coordLocX05 = coordLocX - 0.5;
    const double coordLocY05 = coordLocY - 0.5;
    const double coordLocZ05 = coordLocZ - 0.5;

    const int cellLocX = ngp(coordLocX);
    const int cellLocY = ngp(coordLocY);
    const int cellLocZ = ngp(coordLocZ);
    const int cellLocX05 = ngp(coordLocX05);
    const int cellLocY05 = ngp(coordLocY05);
    const int cellLocZ05 = ngp(coordLocZ05);

    B.x() = fieldB(cellLocX, cellLocY05, cellLocZ05, 0);
    B.y() = fieldB(cellLocX05, cellLocY, cellLocZ05, 1);
    B.z() = fieldB(cellLocX05, cellLocY05, cellLocZ, 2);

    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.5 * dt * q_m * betaI;

    const int blockIndex = sind(cellLocX, cellLocY, cellLocZ);
    auto& currentBlock = LmatX2[blockIndex];

    const int xOffset = cellLocX05 - cellLocX + 1;
    const int yOffset = cellLocY05 - cellLocY + 1;
    const int zOffset = cellLocZ05 - cellLocZ + 1;

    const int indx = BlockDimsNGP::indX(0, yOffset, zOffset);
    const int indy = BlockDimsNGP::indY(xOffset, 0, zOffset);
    const int indz = BlockDimsNGP::indZ(xOffset, yOffset, 0);

    const double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                               {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                               {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    const double common = betaL;
    const int id[3] = {indx, indy, indz};

    for (int c1 = 0; c1 < 3; ++c1) {
        for (int c2 = 0; c2 < 3; ++c2) {
            currentBlock(id[c1], id[c2], c1 * 3 + c2) += common * matB[c1][c2];
        }
    }
}

template void Mesh::update_Lmat2<8>(const Vector3R* coord, int size, const Domain& domain, double charge, double mass,
                                    double mpw, const Field3d& fieldB, const double dt, BlockStack& currentBlock) const;
template void Mesh::update_Lmat2<16>(const Vector3R* coord, int size, const Domain& domain, double charge, double mass,
                                     double mpw, const Field3d& fieldB, const double dt,
                                     BlockStack& currentBlock) const;

template void Mesh::update_Lmat2<32>(const Vector3R* coord, int size, const Domain& domain, double charge, double mass,
                                     double mpw, const Field3d& fieldB, const double dt,
                                     BlockStack& currentBlock) const;

template void Mesh::update_Lmat2<48>(const Vector3R* coord, int size, const Domain& domain, double charge, double mass,
                                     double mpw, const Field3d& fieldB, const double dt,
                                     BlockStack& currentBlock) const;

template void Mesh::update_Lmat2<64>(const Vector3R* coord, int size, const Domain& domain, double charge, double mass,
                                     double mpw, const Field3d& fieldB, const double dt,
                                     BlockStack& currentBlock) const;

// void Mesh::apply_periodic_boundaries(std::vector<IndexMap>& LmatX) {
//     const auto size = Vector3I(xSize, ySize, zSize);
//     constexpr int OVERLAP_SIZE = 3;
//     const int last_indx = size.x() - OVERLAP_SIZE;
//     const int last_indy = size.y() - OVERLAP_SIZE;
//     const int last_indz = size.z() - OVERLAP_SIZE;
//     return;
//     Bounds bounds;
//     if (bounds.isPeriodic(X)) {
// #pragma omp parallel for schedule(dynamic, 32)
//         for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//             auto ix = pos_vind(i, 0);
//             auto iy = pos_vind(i, 1);
//             auto iz = pos_vind(i, 2);
//             auto id = pos_vind(i, 3);

//             if (ix < OVERLAP_SIZE) {
//                 for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it)
//                 {
//                     auto ind2 = it->first;
//                     auto value = it->second;
//                     auto ix1 = pos_vind(ind2, 0);
//                     auto iy1 = pos_vind(ind2, 1);
//                     auto iz1 = pos_vind(ind2, 2);
//                     auto id1 = pos_vind(ind2, 3);
//                     auto indBound = vind(last_indx + ix, iy, iz, id);

//                     if (ix1 < OVERLAP_SIZE) {
//                         auto indBound2 = vind(last_indx + ix1, iy1, iz1,
//                         id1); LmatX[indBound][indBound2] += value;
//                     } else {
//                         LmatX[indBound][ind2] += value;
//                     }
//                 }
//             }
//         }
//     }

//     if (bounds.isPeriodic(Y)) {
// #pragma omp parallel for schedule(dynamic, 32)
//         for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//             auto ix = pos_vind(i, 0);
//             auto iy = pos_vind(i, 1);
//             auto iz = pos_vind(i, 2);
//             auto id = pos_vind(i, 3);

//             if (iy < OVERLAP_SIZE) {
//                 for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it)
//                 {
//                     auto ind2 = it->first;
//                     auto value = it->second;
//                     auto ix1 = pos_vind(ind2, 0);
//                     auto iy1 = pos_vind(ind2, 1);
//                     auto iz1 = pos_vind(ind2, 2);
//                     auto id1 = pos_vind(ind2, 3);
//                     auto indBound = vind(ix, last_indy + iy, iz, id);
//                     if (iy1 < OVERLAP_SIZE) {
//                         auto indBound2 = vind(ix1, last_indy + iy1, iz1,
//                         id1); LmatX[indBound][indBound2] += value;
//                     } else {
//                         LmatX[indBound][ind2] += value;
//                     }
//                 }
//             }
//         }
//     }

//     if (bounds.isPeriodic(Z)) {
// #pragma omp parallel for schedule(dynamic, 32)
//         for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//             auto ix = pos_vind(i, 0);
//             auto iy = pos_vind(i, 1);
//             auto iz = pos_vind(i, 2);
//             auto id = pos_vind(i, 3);

//             if (iz < OVERLAP_SIZE) {
//                 for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it)
//                 {
//                     auto ind2 = it->first;
//                     auto value = it->second;
//                     auto ix1 = pos_vind(ind2, 0);
//                     auto iy1 = pos_vind(ind2, 1);
//                     auto iz1 = pos_vind(ind2, 2);
//                     auto id1 = pos_vind(ind2, 3);
//                     auto indBound = vind(ix, iy, last_indz + iz, id);
//                     if (iz1 < OVERLAP_SIZE) {
//                         auto indBound2 = vind(ix1, iy1, last_indz + iz1,
//                         id1); LmatX[indBound][indBound2] += value;
//                     } else {
//                         LmatX[indBound][ind2] += value;
//                     }
//                 }
//             }
//         }
//     }

//     if (bounds.isPeriodic(X)) {
// #pragma omp parallel for schedule(dynamic, 32)
//         for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//             auto ix = pos_vind(i, 0);
//             auto iy = pos_vind(i, 1);
//             auto iz = pos_vind(i, 2);
//             auto id = pos_vind(i, 3);

//             if (ix > last_indx - 1) {
//                 for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it)
//                 {
//                     auto ind2 = it->first;
//                     auto value = it->second;
//                     auto ix1 = pos_vind(ind2, 0);
//                     auto iy1 = pos_vind(ind2, 1);
//                     auto iz1 = pos_vind(ind2, 2);
//                     auto id1 = pos_vind(ind2, 3);
//                     auto indBound = vind(ix - last_indx, iy, iz, id);

//                     if (ix1 > last_indx - 1) {
//                         auto indBound2 = vind(ix1 - last_indx, iy1, iz1,
//                         id1); LmatX[indBound][indBound2] = value;
//                     } else {
//                         LmatX[indBound][ind2] = value;
//                     }
//                 }
//             }
//         }
//     }

//     if (bounds.isPeriodic(Y)) {
// #pragma omp parallel for schedule(dynamic, 32)
//         for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//             auto ix = pos_vind(i, 0);
//             auto iy = pos_vind(i, 1);
//             auto iz = pos_vind(i, 2);
//             auto id = pos_vind(i, 3);

//             if (iy > last_indy - 1) {
//                 for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it)
//                 {
//                     auto ind2 = it->first;
//                     auto value = it->second;
//                     auto ix1 = pos_vind(ind2, 0);
//                     auto iy1 = pos_vind(ind2, 1);
//                     auto iz1 = pos_vind(ind2, 2);
//                     auto id1 = pos_vind(ind2, 3);
//                     auto indBound = vind(ix, iy - last_indy, iz, id);

//                     if (iy1 > last_indy - 1) {
//                         auto indBound2 = vind(ix1, iy1 - last_indy, iz1,
//                         id1); LmatX[indBound][indBound2] = value;
//                     } else {
//                         LmatX[indBound][ind2] = value;
//                     }
//                 }
//             }
//         }
//     }

//     if (bounds.isPeriodic(Z)) {
// #pragma omp parallel for schedule(dynamic, 32)
//         for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//             auto ix = pos_vind(i, 0);
//             auto iy = pos_vind(i, 1);
//             auto iz = pos_vind(i, 2);
//             auto id = pos_vind(i, 3);

//             if (iz > last_indz - 1) {
//                 for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it)
//                 {
//                     auto ind2 = it->first;
//                     auto value = it->second;
//                     auto ix1 = pos_vind(ind2, 0);
//                     auto iy1 = pos_vind(ind2, 1);
//                     auto iz1 = pos_vind(ind2, 2);
//                     auto id1 = pos_vind(ind2, 3);
//                     auto indBound = vind(ix, iy, iz - last_indz, id);

//                     if (iz1 > last_indz - 1) {
//                         auto indBound2 = vind(ix1, iy1, iz1 - last_indz,
//                         id1); LmatX[indBound][indBound2] = value;
//                     } else {
//                         LmatX[indBound][ind2] = value;
//                     }
//                 }
//             }
//         }
//     }
// }
