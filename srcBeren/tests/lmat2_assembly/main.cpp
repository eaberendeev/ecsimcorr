// Lmat2 matrix assembly equivalence test.
//
// Compares 5 independent assembly paths of the implicit ECSIM matrix Lmat2 on
// many particle scenarios:
//   P1 Reference  : fill_matrixL2 + stencil_Lmat2_Reference
//   P2 Old optim. : fill_matrixL2 + stencil_Lmat2
//   P3 New V2     : fill_matrixL2_Optimized + stencil_Lmat2_Optimized_V2
//   P4 Dict oracle A: independent re-read of mesh.LmatX2 blocks -> std::map
//   P5 Dict oracle B: direct per-particle rebuild from coordinates/physics
//                     (independent of LmatX2 blocks AND BlockDims indexing)
//
// All paths must agree in structure (|v| <= 1e-14 treated as zero) and in
// values (rel. tol 1e-9, abs. tol 1e-14 for near-zeros).

#include <omp.h>

#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "Mesh.h"
#include "Particle.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "bmatrix.h"
#include "nlohmann/json.hpp"
#include "row_block.h"
#include "types.h"

// =============================================================================
// tiny test framework (same style as srcBeren/tests/unit/domain.cpp)
// =============================================================================
namespace test {
static int passed = 0;
static int failed = 0;

void assert_true(bool condition, const std::string& test_name) {
    if (condition) {
        passed++;
    } else {
        std::cout << "[FAIL] " << test_name << std::endl;
        failed++;
    }
}

void print_summary() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

int get_failed() {
    return failed;
}
}   // namespace test

// =============================================================================
// configuration constants
// =============================================================================
constexpr int NCELL = 6;            // cells per dimension
constexpr double CS = 0.1;          // cell size
constexpr double DT = 0.1;          // time step
constexpr double DICT_TOL = 1e-16;  // block-entry tolerance (same as blockToTriplets TOL)
constexpr double ZERO_TOL = 1e-14;  // structure filter: |v| <= ZERO_TOL is "zero"
constexpr double REL_TOL = 1e-9;    // relative value tolerance
constexpr double ABS_TOL = 1e-14;   // absolute value tolerance
constexpr int MAX_MISMATCH_PRINT = 10;
constexpr double SANITY_MIN_ABS_VALUE = 1e-6;  // sanity: matrix values must be non-trivial

// =============================================================================
// sparse-matrix comparison
// =============================================================================
struct Entry {
    int row;
    int col;
    double val;
};

static bool operator<(const Entry& a, const Entry& b) {
    if (a.row != b.row)
        return a.row < b.row;
    return a.col < b.col;
}

// Extract (row, col, val) triplets skipping entries with |val| <= tol (treated
// as structurally zero). Sorted by (row, col).
static std::vector<Entry> extractEntries(const Operator& mat, double tol) {
    std::vector<Entry> res;
    for (int row = 0; row < mat.outerSize(); ++row) {
        for (Operator::InnerIterator it(mat, row); it; ++it) {
            if (std::abs(it.value()) > tol) {
                res.push_back({row, static_cast<int>(it.col()), it.value()});
            }
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}

static bool valuesClose(double a, double b) {
    const double diff = std::fabs(a - b);
    return diff <= ABS_TOL + REL_TOL * std::max(std::fabs(a), std::fabs(b));
}

// Compare two CSR matrices: same structure after tolerance filtering, values
// within tolerance. Prints the first MAX_MISMATCH_PRINT mismatches.
static bool compareCSR(const Operator& matA, const Operator& matB, const std::string& label) {
    if (matA.rows() != matB.rows() || matA.cols() != matB.cols()) {
        std::cout << "  " << label << ": dimension mismatch " << matA.rows() << "x" << matA.cols() << " vs "
                  << matB.rows() << "x" << matB.cols() << std::endl;
        return false;
    }

    const std::vector<Entry> a = extractEntries(matA, ZERO_TOL);
    const std::vector<Entry> b = extractEntries(matB, ZERO_TOL);

    if (a.size() != b.size()) {
        std::cout << "  " << label << ": nnz mismatch after tol filter: A=" << a.size() << " B=" << b.size()
                  << std::endl;
    }

    bool ok = (a.size() == b.size());
    int printed = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i].row == b[j].row && a[i].col == b[j].col) {
            if (!valuesClose(a[i].val, b[j].val)) {
                if (printed < MAX_MISMATCH_PRINT) {
                    std::cout << "  " << label << ": value mismatch at (row " << a[i].row << ", col " << a[i].col
                              << "): A=" << a[i].val << " B=" << b[j].val << std::endl;
                    printed++;
                }
                ok = false;
            }
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            if (printed < MAX_MISMATCH_PRINT) {
                std::cout << "  " << label << ": entry only in A at (row " << a[i].row << ", col " << a[i].col
                          << "): A=" << a[i].val << std::endl;
                printed++;
            }
            ok = false;
            ++i;
        } else {
            if (printed < MAX_MISMATCH_PRINT) {
                std::cout << "  " << label << ": entry only in B at (row " << b[j].row << ", col " << b[j].col
                          << "): B=" << b[j].val << std::endl;
                printed++;
            }
            ok = false;
            ++j;
        }
    }
    while (i < a.size()) {
        if (printed < MAX_MISMATCH_PRINT) {
            std::cout << "  " << label << ": entry only in A at (row " << a[i].row << ", col " << a[i].col
                      << "): A=" << a[i].val << std::endl;
            printed++;
        }
        ok = false;
        ++i;
    }
    while (j < b.size()) {
        if (printed < MAX_MISMATCH_PRINT) {
            std::cout << "  " << label << ": entry only in B at (row " << b[j].row << ", col " << b[j].col
                      << "): B=" << b[j].val << std::endl;
            printed++;
        }
        ok = false;
        ++j;
    }

    if (printed >= MAX_MISMATCH_PRINT && (a.size() != b.size() || ok == false)) {
        std::cout << "  " << label << ": ... further mismatches not shown" << std::endl;
    }
    if (ok) {
        std::cout << "  " << label << ": OK (nnz = " << a.size() << ")" << std::endl;
    }
    return ok;
}

// =============================================================================
// dict (std::map keyed by global row) vs CSR comparison
// =============================================================================
using Dict = std::map<int, std::map<int, double>>;

static int dictSize(const Dict& dict) {
    int n = 0;
    for (const auto& kv : dict) {
        n += static_cast<int>(kv.second.size());
    }
    return n;
}

// Materialize the dict oracle as an Eigen CSR matrix via setFromTriplets —
// the same trusted accumulation path the P1 reference uses. Duplicate
// (row, col) triplets are summed by Eigen, matching the dict accumulation.
static Operator dictToCSR(const Dict& dict, int n) {
    std::vector<Trip> trips;
    trips.reserve(dictSize(dict));
    for (const auto& [row, cols] : dict) {
        for (const auto& [col, val] : cols) {
            if (std::abs(val) > ZERO_TOL) {
                trips.emplace_back(row, col, val);
            }
        }
    }
    Operator mat(n, n);
    mat.setFromTriplets(trips.begin(), trips.end());
    return mat;
}


// =============================================================================
// P4: independent dictionary oracle from mesh.LmatX2 blocks
//
// Re-implements the neighbor-cell / indexer mapping used inside
// blockToTriplets<XIndexer, YIndexer, d> in operators.cpp (~line 178)
// WITHOUT calling it: reads block values directly via
// calculateIndex(i, j, d) = j + 12 * (12 * d + i) (bmatrix.h) and accumulates
// into a std::map<int, std::map<int,double>>.
// =============================================================================
static Dict buildDictFromBlocks(const Mesh& mesh, const Domain& domain) {
    Dict dict;
    const int Ny = domain.size().y();
    const int Nz = domain.size().z();

    // vind with the same formula as blockToTriplets (no asserts, allows ghost
    // index 0 even if a node sits outside the physical box).
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };

    // Indexer tables replicating XIndexer/YIndexer/ZIndexer from bmatrix.h:
    //   X: 3x2x2, dir 0, offset (-1, 0, 0)
    //   Y: 2x3x2, dir 1, offset (0, -1, 0)
    //   Z: 2x2x3, dir 2, offset (0, 0, -1)
    struct Idx {
        int dir, offx, offy, offz;
        int sx, sy, sz;
    };
    static constexpr Idx T[3] = {
        {0, -1, 0, 0, 3, 2, 2},
        {1, 0, -1, 0, 2, 3, 2},
        {2, 0, 0, -1, 2, 2, 3},
    };
    // Indexer::calculate(x, y, z) = x * (sy * sz) + y * sz + z
    auto calc = [](int c, int x, int y, int z) -> int {
        const Idx& t = T[c];
        return x * (t.sy * t.sz) + y * t.sz + z;
    };

    constexpr int BORDER = 1;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;

    for (int i = BORDER; i < max_i; ++i) {
        for (int j = BORDER; j < max_j; ++j) {
            for (int k = BORDER; k < max_k; ++k) {
                const int sind = mesh.sind(i, j, k);
                if (!mesh.LmatX2.non_zeros[sind]) {
                    continue;
                }
                const Block& block = mesh.LmatX2[sind];
                for (int c1 = 0; c1 < 3; ++c1) {
                    for (int c2 = 0; c2 < 3; ++c2) {
                        const int DIR = c1 * 3 + c2;
                        const Idx& R = T[c1];
                        const Idx& C = T[c2];
                        for (int x1 = 0; x1 < R.sx; ++x1) {
                            for (int y1 = 0; y1 < R.sy; ++y1) {
                                for (int z1 = 0; z1 < R.sz; ++z1) {
                                    const int row =
                                        vind(i + x1 + R.offx, j + y1 + R.offy, k + z1 + R.offz, R.dir);
                                    const int rowLocal = calc(c1, x1, y1, z1);
                                    for (int x2 = 0; x2 < C.sx; ++x2) {
                                        for (int y2 = 0; y2 < C.sy; ++y2) {
                                            for (int z2 = 0; z2 < C.sz; ++z2) {
                                                const int colLocal = calc(c2, x2, y2, z2);
                                                // calculateIndex(rowLocal, colLocal, DIR)
                                                const double val =
                                                    block.values[colLocal + 12 * (12 * DIR + rowLocal)];
                                                if (std::abs(val) > DICT_TOL) {
                                                    const int col =
                                                        vind(i + x2 + C.offx, j + y2 + C.offy, k + z2 + C.offz,
                                                             C.dir);
                                                    dict[row][col] += val;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return dict;
}

// =============================================================================
// P5: direct-from-particles dictionary oracle (oracle B)
//
// Rebuilds the Lmat2 matrix from first principles, replicating the SEMANTICS
// of Mesh::update_Lmat2 (srcBeren/fields/mesh/Mesh.cpp:81-186) but with a
// completely different structure: instead of accumulating particle
// contributions into LmatX2 blocks addressed via BlockDims::indX/indY/indZ,
// every particle directly emits (global row, global col, value) triplets that
// are summed into a std::map.  This verifies the whole chain: per-cell block
// accumulation physics + index mapping + final CSR assembly.
//
// Derivation of the block-local index -> global vind mapping, INDEPENDENT of
// the BlockDims / XIndexer tables (verified against them, see below):
//
// update_Lmat2 writes, for the component pair (c1, c2) and a local shape
// sample (i,j,k)x(i1,j1,k1):
//   block(rowIndex, colIndex, c1*3+c2) += betaL * s1[c1] * s2[c2] * matB[c1][c2]
// with (xOffset = cellLocX05 - cellLocX + 1, likewise y/z):
//   rowIndex = indX(xOffset+i, j, k)   for c1=0 (X)   = 4*(xOffset+i) + 2*j + k
//              indY(i, yOffset+j, k)   for c1=1 (Y)   = 6*i + 2*(yOffset+j) + k
//              indZ(i, j, zOffset+k)   for c1=2 (Z)   = 6*i + 3*j + (zOffset+k)
// and analogously for colIndex with (i1,j1,k1).
//
// blockToTriplets (operators.cpp:177-199) converts the local XIndexer/
// YIndexer/ZIndexer node (x1,y1,z1) of the block of cell (cx,cy,cz) to the
// global E* grid node via (local index equals indX/indY/indZ respectively):
//   X: node=(cx-1+x1, cy+y1, cz+z1), dir 0, local index 4*x1 + 2*y1 + z1
//   Y: node=(cx+x1,   cy-1+y1, cz+z1), dir 1, local index 6*x1 + 2*y1 + z1
//   Z: node=(cx+x1,   cy+y1,   cz-1+z1), dir 2, local index 6*x1 + 3*y1 + z1
// (the XIndexer 3x2x2 / YIndexer 2x3x2 / ZIndexer 2x2x3 sizes match the
// BlockDims::X/Y/Z SIZE_* constants, and XIndexer::calculate == indX, etc.)
//
// Equating the two local-index expressions:
//   X: 4*(xOffset+i)+2j+k = 4*x1+2*y1+z1  =>  x1=xOffset+i, y1=j, z1=k
//   Y: 6*i+2*(yOffset+j)+k = 6*x1+2*y1+z1 =>  x1=i, y1=yOffset+j, z1=k
//   Z: 6*i+3*j+(zOffset+k) = 6*x1+3*y1+z1 =>  x1=i, y1=j, z1=zOffset+k
// and substituting xOffset = cellLocX05 - cellLocX + 1 collapses the global
// node's x index to cellLocX - 1 + x1 = cellLocX05 + i (likewise y/z) giving
// the DIRECT formula used below (independent of both BlockDims and Indexer):
//   c1=0 (Ex): global node = (cellLocX05+i, cellLocY+j,   cellLocZ+k)
//   c1=1 (Ey): global node = (cellLocX+i,   cellLocY05+j, cellLocZ+k)
//   c1=2 (Ez): global node = (cellLocX+i,   cellLocY+j,   cellLocZ05+k)
//   global row = vind(node, c1);  global col = vind(node', c2) with (i1,j1,k1).
//
// s1/s2 weights are the component-row tuples computed in update_Lmat2:
//   s1[0]=sx05[i]*sy[j]*sz[k], s1[1]=sx[i]*sy05[j]*sz[k], s1[2]=sx[i]*sy[j]*sz05[k]
// (s2 identical in (i1,j1,k1)).
// =============================================================================
static Dict buildDictFromParticles(const ParticlesArray& sp, const Field3d& fieldB, const Domain& domain,
                                   double dt) {
    Dict dict;
    const int Ny = domain.size().y();
    const int Nz = domain.size().z();

    // vind with the same formula as blockToTriplets (no asserts, allows ghost
    // index 0 even if a node sits outside the physical box).
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };

    constexpr int SMAX = 2;   // SHAPE_SIZE
    const double q_m = sp.charge / sp.mass();

    // Iterate particlesData the same way fill_matrixL_impl_linear2 does
    // (ParticlesCore.cpp:72-83): per cell, per particle, its own coord.
    for (int pk = 0; pk < sp.particlesData.capacity(); ++pk) {
        const std::vector<Particle>& cell = sp.particlesData(pk);
        if (cell.empty()) {
            continue;
        }
        for (const Particle& p : cell) {
            const Vector3R& coord = p.coord;   // Lmat2 has NO velocity dependence

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

            double sx[SMAX], sy[SMAX], sz[SMAX], sx05[SMAX], sy05[SMAX], sz05[SMAX];
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

            // B field at the particle position, Yee staggered (update_Lmat2
            // Mesh.cpp:117-134).
            Vector3R B(0.);
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
                        B.x() += wx * fieldB(indx, indy05, indz05, 0);
                        B.y() += wy * fieldB(indx05, indy, indz05, 1);
                        B.z() += wz * fieldB(indx05, indy05, indz, 2);
                    }
                }
            }

            const Vector3R b = 0.5 * dt * q_m * B;
            const double betaI = sp.mpw() * sp.charge / (1.0 + b.squared());
            const double betaL = 0.25 * dt * dt * q_m * betaI;

            const double matB[3][3] = {{1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
                                       {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
                                       {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

            // Direct global-node vind + component linear weight for the local
            // shape sample (i,j,k) of component c (derivation above).
            auto nodeVind = [&](int c, int i, int j, int k) -> int {
                switch (c) {
                    case 0:
                        return vind(cellLocX05 + i, cellLocY + j, cellLocZ + k, 0);
                    case 1:
                        return vind(cellLocX + i, cellLocY05 + j, cellLocZ + k, 1);
                    default:
                        return vind(cellLocX + i, cellLocY + j, cellLocZ05 + k, 2);
                }
            };
            auto shape = [&](int c, int i, int j, int k) -> double {
                switch (c) {
                    case 0:
                        return sx05[i] * sy[j] * sz[k];
                    case 1:
                        return sx[i] * sy05[j] * sz[k];
                    default:
                        return sx[i] * sy[j] * sz05[k];
                }
            };

            for (int c1 = 0; c1 < 3; ++c1) {
                for (int c2 = 0; c2 < 3; ++c2) {
                    const double core = betaL * matB[c1][c2];
                    for (int i = 0; i < SMAX; ++i) {
                        for (int j = 0; j < SMAX; ++j) {
                            for (int k = 0; k < SMAX; ++k) {
                                const int row = nodeVind(c1, i, j, k);
                                const double s1 = shape(c1, i, j, k);
                                for (int i1 = 0; i1 < SMAX; ++i1) {
                                    for (int j1 = 0; j1 < SMAX; ++j1) {
                                        for (int k1 = 0; k1 < SMAX; ++k1) {
                                            const double val = core * s1 * shape(c2, i1, j1, k1);
                                            if (std::abs(val) > DICT_TOL) {
                                                dict[row][nodeVind(c2, i1, j1, k1)] += val;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return dict;
}

// =============================================================================
// assembly helpers
// =============================================================================
struct AssemblyResult {
    Operator matRef;   // P1
    Operator matOld;   // P2
    Operator matV2;    // P3
    Dict dict;         // P4
    Dict dictB;        // P5
    long long nnzRef = 0;
};

static nlohmann::json makeSpeciesConfig() {
    nlohmann::json cfg;
    cfg["Name"] = "lmat2_test_species";
    cfg["Charge"] = -1.0;
    cfg["Density"] = 1.0;
    cfg["Mass"] = 1.0;
    cfg["NumPartPerCell"] = 10;
    return cfg;
}

// Prepare (zero) the LmatX2 blocks with the same semantics as
// SimulationEcsim::prepare_block_matrix (simulation_ecsim.cpp:164).
static void prepareLmatX2(Mesh& mesh, const ParticlesArray& sp) {
    Array3D<int> countInCell(sp.particlesData.size());
    countInCell.setZero();
    for (int i = 0; i < sp.particlesData.capacity(); ++i) {
        countInCell(i) += static_cast<int>(sp.particlesData(i).size());
    }
    mesh.LmatX2.prepare(countInCell);
}

static void fillUniformFieldB(Field3d& fieldB) {
    const Vector3I size = fieldB.sizes();
    for (int i = 0; i < size.x(); ++i) {
        for (int j = 0; j < size.y(); ++j) {
            for (int k = 0; k < size.z(); ++k) {
                fieldB(i, j, k, 0) = 0.1 + 0.001 * (i + j + k);   // mildly non-uniform
                fieldB(i, j, k, 1) = 0.2;
                fieldB(i, j, k, 2) = 0.3;
            }
        }
    }
}

// Runs all 4 assembly paths on a fresh/clean mesh state and compares them.
// Requires LmatX2 already prepared (by caller) and rowBlocksGlobal resized(0).
static void assembleAndCompare(const std::string& scenario, const std::string& stage, Mesh& mesh,
                               ParticlesArray& sp, const Field3d& fieldB, const Domain& domain,
                               std::vector<std::vector<RowBlock<36>>>& rowBlocksGlobal, AssemblyResult& out) {
    const int n = domain.total_size() * 3;

    // -- shared LmatX2 fill for P1 / P2 / P4 --------------------------------
    sp.fill_matrixL2(mesh, fieldB, domain, DT, ShapeType::Linear);

    // P1: reference
    out.matRef.resize(n, n);
    mesh.stencil_Lmat2_Reference(out.matRef, domain);

    // P2: old optimized
    out.matOld.resize(n, n);
    mesh.stencil_Lmat2(out.matOld, domain, mesh.workspacePtr);

    // P4: independent dict oracle on the same LmatX2 blocks
    out.dict = buildDictFromBlocks(mesh, domain);

    // P5: direct per-particle dict oracle (independent of LmatX2 and
    // BlockDims indexing).
    out.dictB = buildDictFromParticles(sp, fieldB, domain, DT);

    // P3: new V2 (production path, rowBlocksGlobal must be reset before use)
    for (int i = 0; i < static_cast<int>(rowBlocksGlobal.size()); ++i) {
        rowBlocksGlobal[i].resize(0);
        rowBlocksGlobal[i].reserve(1024 * 1024 * 4);
    }
    sp.fill_matrixL2_Optimized(mesh, fieldB, domain, DT, ShapeType::Linear, rowBlocksGlobal);
    out.matV2.resize(n, n);
    mesh.stencil_Lmat2_Optimized_V2(out.matV2, domain, rowBlocksGlobal, mesh.workspacePtr);

    out.nnzRef = extractEntries(out.matRef, ZERO_TOL).size();

    std::cout << "  [" << stage << "] particles=" << sp.get_total_num_of_particles() << " nnz(P1)=" << out.nnzRef
              << " nnz(P5)=" << dictSize(out.dictB) << std::endl;

    // -- sanity: inputs and results must be non-trivial ---------------------
    // The matrix value is betaL * s1*s2 * matB with
    //   betaL = 0.25*dt^2*(q/m)*(mpw*q/(1+b^2)); it depends on coord, B,
    //   q/m, mpw, dt (particle velocity does NOT enter Lmat2). If any input
    //   silently parsed to zero, all 4 paths would agree on an all-zero
    //   matrix and the comparisons would pass vacuously. Guard against that.
    test::assert_true(sp.charge != 0.0, scenario + " [" + stage + "] sanity: charge != 0");
    test::assert_true(sp.mass() != 0.0, scenario + " [" + stage + "] sanity: mass != 0");
    test::assert_true(sp.mpw() > 0.0, scenario + " [" + stage + "] sanity: mpw > 0");
    double maxFieldB = 0.0;
    for (int i = 0; i < fieldB.sizes().x(); ++i)
        for (int j = 0; j < fieldB.sizes().y(); ++j)
            for (int k = 0; k < fieldB.sizes().z(); ++k)
                for (int d = 0; d < 3; ++d) maxFieldB = std::max(maxFieldB, std::fabs(fieldB(i, j, k, d)));
    test::assert_true(maxFieldB > 0.0, scenario + " [" + stage + "] sanity: fieldB != 0");
    // non-empty scenarios must produce a non-empty matrix with non-trivial values
    if (sp.get_total_num_of_particles() > 0) {
        test::assert_true(out.nnzRef > 0, scenario + " [" + stage + "] sanity: nnz(P1) > 0");
        double maxAbs = 0.0;
        for (int row = 0; row < out.matRef.outerSize(); ++row)
            for (Operator::InnerIterator it(out.matRef, row); it; ++it)
                maxAbs = std::max(maxAbs, std::fabs(it.value()));
        test::assert_true(maxAbs > SANITY_MIN_ABS_VALUE,
                          scenario + " [" + stage + "] sanity: max|Lmat2| > " + std::to_string(SANITY_MIN_ABS_VALUE) +
                              " (got " + std::to_string(maxAbs) + ")");
    }

    const std::string p = scenario + " [" + stage + "] P1==P2";
    const std::string p2 = scenario + " [" + stage + "] P1==P3";
    const std::string p3 = scenario + " [" + stage + "] P2==P3";
    const std::string d4 = scenario + " [" + stage + "] P1==P4(dict)";
    const std::string d5 = scenario + " [" + stage + "] P2==P4(dict)";
    const std::string d6 = scenario + " [" + stage + "] P3==P4(dict)";
    const std::string d7 = scenario + " [" + stage + "] P1==P5(dictB)";
    const std::string d8 = scenario + " [" + stage + "] P2==P5(dictB)";
    const std::string d9 = scenario + " [" + stage + "] P3==P5(dictB)";

    test::assert_true(compareCSR(out.matRef, out.matOld, "P1 vs P2"), p);
    test::assert_true(compareCSR(out.matRef, out.matV2, "P1 vs P3"), p2);
    test::assert_true(compareCSR(out.matOld, out.matV2, "P2 vs P3"), p3);
    // P4 dict is materialized as CSR via setFromTriplets (the same trusted
    // Eigen accumulation path as the P1 reference) and compared with the
    // single unified comparator.
    Operator matDict = dictToCSR(out.dict, n);
    test::assert_true(compareCSR(out.matRef, matDict, "P1 vs P4(dict)"), d4);
    test::assert_true(compareCSR(out.matOld, matDict, "P2 vs P4(dict)"), d5);
    test::assert_true(compareCSR(out.matV2, matDict, "P3 vs P4(dict)"), d6);
    // P5 dictB likewise, against all three production CSR paths.
    Operator matDictB = dictToCSR(out.dictB, n);
    test::assert_true(compareCSR(out.matRef, matDictB, "P1 vs P5(dictB)"), d7);
    test::assert_true(compareCSR(out.matOld, matDictB, "P2 vs P5(dictB)"), d8);
    test::assert_true(compareCSR(out.matV2, matDictB, "P3 vs P5(dictB)"), d9);
}

// =============================================================================
// scenario helpers
// =============================================================================
static Domain makeDomain() {
    Domain domain;
    domain.init(Vector3I(NCELL, NCELL, NCELL), Vector3R(CS, CS, CS));
    return domain;
}

static double uniform(std::mt19937& gen, double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(gen);
}

// Random displacement used by the reassembly scenario.
static Vector3R randomDisplacement(std::mt19937& gen) {
    return Vector3R(uniform(gen, -0.03, 0.03), uniform(gen, -0.03, 0.03), uniform(gen, -0.03, 0.03));
}

static void runScenario(const std::string& scenario, void (*populate)(ParticlesArray&, std::mt19937&, bool)) {
    Domain domain = makeDomain();
    BoundaryConditionHandler bc_handler;   // empty handler: no BC modifications
    Mesh mesh;
    mesh.init(domain, DT, bc_handler);

    Field3d fieldB(domain.size(), 3);
    fieldB.setZero();
    fillUniformFieldB(fieldB);

    ParticlesArray sp(makeSpeciesConfig(), domain);
    std::mt19937 gen(20260601);
    populate(sp, gen, false);
    sp.update_cells(domain);   // groups particles into cells (required by P3)

    prepareLmatX2(mesh, sp);

    std::cout << "--- scenario: " << scenario << " ---" << std::endl;
    std::vector<std::vector<RowBlock<36>>> rowBlocksGlobal(omp_get_max_threads());

    AssemblyResult out;
    assembleAndCompare(scenario, "assemble-1", mesh, sp, fieldB, domain, rowBlocksGlobal, out);

    if (scenario == "reassembly") {
        // move particles a little, regroup, re-prepare blocks and re-assemble
        // into the SAME Operator/buffers to catch stale state.
        std::mt19937 gen(20260707);
        for (int i = 0; i < sp.particlesData.capacity(); ++i) {
            for (auto& p : sp.particlesData(i)) {
                Vector3R d = randomDisplacement(gen);
                p.coord += d;
                p.coord.x() = std::min(std::max(p.coord.x(), 0.002), 0.598);
                p.coord.y() = std::min(std::max(p.coord.y(), 0.002), 0.598);
                p.coord.z() = std::min(std::max(p.coord.z(), 0.002), 0.598);
            }
        }
        sp.update_cells(domain);
        prepareLmatX2(mesh, sp);
        // old operators and rowBlocksGlobal buffers are reused as-is
        assembleAndCompare(scenario, "assemble-2", mesh, sp, fieldB, domain, rowBlocksGlobal, out);
    }
}

// ---------------------------- populate functions ----------------------------
static void popEmpty(ParticlesArray& sp, std::mt19937&, bool) {
    (void) sp;
}

static void popSingleCenter(ParticlesArray& sp, std::mt19937&, bool) {
    // one particle exactly at the center of an interior cell
    sp.add_particle(Particle(Vector3R(0.25, 0.25, 0.25), Vector3R(0.0, 0.0, 0.0)));
}

static void popSingleEdge(ParticlesArray& sp, std::mt19937&, bool) {
    // one particle exactly on cell-boundary coordinates (x = 0.1 = one cell size)
    sp.add_particle(Particle(Vector3R(0.1, 0.2, 0.3), Vector3R(0.0, 0.0, 0.0)));
}

static void popAllInOneCell(ParticlesArray& sp, std::mt19937& gen, bool) {
    // 10 particles with tiny jitter inside one interior cell (center 0.25^3)
    for (int n = 0; n < 10; ++n) {
        const double dx = uniform(gen, -0.005, 0.005);
        const double dy = uniform(gen, -0.005, 0.005);
        const double dz = uniform(gen, -0.005, 0.005);
        sp.add_particle(Particle(Vector3R(0.25 + dx, 0.25 + dy, 0.25 + dz),
                                 Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
    }
}

static void popFullFill(ParticlesArray& sp, std::mt19937& gen, bool) {
    // 10 random particles in EVERY cell (full matrix coverage)
    for (int ci = 0; ci < NCELL; ++ci) {
        for (int cj = 0; cj < NCELL; ++cj) {
            for (int ck = 0; ck < NCELL; ++ck) {
                const double cx = (ci + 0.5) * CS;
                const double cy = (cj + 0.5) * CS;
                const double cz = (ck + 0.5) * CS;
                for (int n = 0; n < 10; ++n) {
                    const double dx = uniform(gen, -0.49 * CS, 0.49 * CS);
                    const double dy = uniform(gen, -0.49 * CS, 0.49 * CS);
                    const double dz = uniform(gen, -0.49 * CS, 0.49 * CS);
                    sp.add_particle(Particle(Vector3R(cx + dx, cy + dy, cz + dz),
                                             Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1),
                                                      uniform(gen, -0.1, 0.1))));
                }
            }
        }
    }
}

static void popRandomScatter(ParticlesArray& sp, std::mt19937& gen, bool) {
    // ~200 particles uniformly random over the domain [0, 0.6)^3
    const int count = 200;
    for (int n = 0; n < count; ++n) {
        sp.add_particle(Particle(Vector3R(uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6)),
                                 Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
    }
}

static void popBoundaries(ParticlesArray& sp, std::mt19937& gen, bool) {
    // particles clustered near domain faces/corners and ghost-adjacent cells.
    // Physical cells map to grid indices {1..6}; index 1 and 6 touch ghosts.
    const double nearFace[] = {0.001, 0.599};
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
            for (int c = 0; c < 2; ++c) {
                const double cx = nearFace[a];
                const double cy = nearFace[b];
                const double cz = nearFace[c];
                for (int n = 0; n < 8; ++n) {
                    const double dx = uniform(gen, -0.02, 0.02);
                    const double dy = uniform(gen, -0.02, 0.02);
                    const double dz = uniform(gen, -0.02, 0.02);
                    sp.add_particle(Particle(Vector3R(std::min(std::max(cx + dx, 0.0), 0.599),
                                                      std::min(std::max(cy + dy, 0.0), 0.599),
                                                      std::min(std::max(cz + dz, 0.0), 0.599)),
                                             Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1),
                                                      uniform(gen, -0.1, 0.1))));
                }
            }
        }
    }
    // a handful exactly on the box corner
    sp.add_particle(Particle(Vector3R(0.0, 0.0, 0.0), Vector3R(0.0, 0.0, 0.0)));
    sp.add_particle(Particle(Vector3R(0.599, 0.599, 0.599), Vector3R(0.0, 0.0, 0.0)));
    sp.add_particle(Particle(Vector3R(0.0, 0.599, 0.0), Vector3R(0.0, 0.0, 0.0)));
    sp.add_particle(Particle(Vector3R(0.599, 0.0, 0.599), Vector3R(0.0, 0.0, 0.0)));
}

static void popReassembly(ParticlesArray& sp, std::mt19937& gen, bool) {
    // like random_scatter; the second assembly (after displacement) is handled
    // inside runScenario for the "reassembly" tag.
    for (int n = 0; n < 120; ++n) {
        sp.add_particle(Particle(Vector3R(uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6)),
                                 Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
    }
}

// =============================================================================
// main
// =============================================================================
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Running Lmat2 Assembly Equivalence Tests" << std::endl;
    std::cout << "OMP_NUM_THREADS = " << omp_get_max_threads() << std::endl;
    std::cout << "========================================\n" << std::endl;

    runScenario("empty", popEmpty);
    runScenario("single_center", popSingleCenter);
    runScenario("single_edge", popSingleEdge);
    runScenario("all_in_one_cell", popAllInOneCell);
    runScenario("full_fill", popFullFill);
    runScenario("random_scatter", popRandomScatter);
    runScenario("boundaries", popBoundaries);
    runScenario("reassembly", popReassembly);

    test::print_summary();
    return test::get_failed();
}