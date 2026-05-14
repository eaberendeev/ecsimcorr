#include "Mesh.h"

#include "Shape.h"
#include "World.h"
#include "interpolation.h"
#include "operators.h"
#include "solverSLE.h"
#include "timer.h"

void Mesh::init(const Domain& domain, double dt) {
    Lmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Lmat2.resize(domain.total_size() * 3, domain.total_size() * 3);
    Mmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Imat.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlE.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlB.resize(domain.total_size() * 3, domain.total_size() * 3);
    IMmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    chargeDensityOld.resize(domain.size(), 1);
    chargeDensity.resize(domain.size(), 1);
    divE.resize(domain.total_size(), domain.total_size() * 3);
    // TODO: move sind, converter func to BlockMatrix, resize with 3dim
    LmatX2.resize(domain.total_size());
    LmatX_NGP.resize(domain.total_size());
    // LmatX2.reserve();

    xCellSize = domain.cell_size().x();
    yCellSize = domain.cell_size().y();
    zCellSize = domain.cell_size().z();
    xSize = domain.size().x();
    ySize = domain.size().y();
    zSize = domain.size().z();

    stencil_Imat(Imat, domain);
    stencil_curlE(curlE, domain);
    stencil_curlB(curlB, domain);

    stencil_divE(divE, domain);

    Mmat = -0.25 * dt * dt * curlB * curlE;
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

// Solve Ax=b for find fieldE.
// (E_{n+1} - E_n) / dt = -J_{n+1/2} + rot(B_{n+1/2}) B_{n+1/2} =
// (B_n + B_{n+1})/2 (B_{n+1}- B_n) / dt = - rot(E_{n+1/2})
// M = -0.25 * dt * dt * rot_opB * rot_opE;
// E_{n+1} = dt*E_n + M*(E_{n+1}+E_n) - dt*J_{n+1/2}+ rot(B_n)

// E_{n+1} - fieldEnew (out)
// E_n - fieldE (in)
// B_n - fieldB (in)
// J_{n+1/2} - fieldJ (in)
void Mesh::impicit_find_fieldE(Field3d& Enew, const Field3d& E, const Field3d& B, const Field3d& J, const double dt) {
    Field rhs = E.data() - dt * J.data() + dt * curlB * B.data() + Mmat * E.data();
    Operator A = Imat - Mmat;
    // TODO: use it for Field3d
    // solve_linear_system<BicgstabSolver<Field>>(
    //     A, rhs, Enew.data(), E.data());

    std::cout << "Solver impicit_find_fieldE error = " << (A * Enew.data() - rhs).norm() << "\n";
}

double Mesh::calculate_residual(const Field3d& Enew, const Field3d& E, const Field3d& B, const Field3d& J,
                                const double dt) {
    Field rhs = E.data() - dt * J.data() + dt * curlB * B.data() + Mmat * E.data();
    Operator A = Imat - Mmat;

    return (A * Enew.data() - rhs).norm();
}

void Mesh::fdtd_explicit(Field3d& E, Field3d& B, const Field3d& J, const double dt) {
    E.data() += 0.5 * dt * curlB * B.data() - 0.5 * dt * J.data();
    B.data() -= 0.5 * dt * curlE * E.data();
}

void Mesh::computeB(const Field3d& fieldE, const Field3d& fieldEn, Field3d& fieldB, double dt) {
    fieldB.data() -= 0.5 * dt * curlE * (fieldE.data() + fieldEn.data());
}

void Mesh::compute_fieldB(Field3d& Bn, const Field3d& B, const Field3d& E, const Field3d& En, double dt) {
    Bn.data() = B.data() - 0.5 * dt * curlE * (E.data() + En.data());
}

void Mesh::update_Lmat2(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                        const Field3d& fieldB, const double dt) {
    RECORD_TIMER;

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
    auto& currentBlock = LmatX2[blockIndex];

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
                const int idx1[3] = {indX(xOffset + i, j, k), indY(i, yOffset + j, k), indZ(i, j, zOffset + k)};

                for (int i1 = 0; i1 < SMAX; ++i1) {
                    for (int j1 = 0; j1 < SMAX; ++j1) {
                        for (int k1 = 0; k1 < SMAX; ++k1) {
                            const double s2[3] = {sx05[i1] * sy[j1] * sz[k1], sx[i1] * sy05[j1] * sz[k1],
                                                  sx[i1] * sy[j1] * sz05[k1]};
                            const int idx2[3] = {indX(xOffset + i1, j1, k1), indY(i1, yOffset + j1, k1),
                                                 indZ(i1, j1, zOffset + k1)};

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

void Mesh::update_Lmat2_NGP(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                            const Field3d& fieldB, const double dt) {
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
