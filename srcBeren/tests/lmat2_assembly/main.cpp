// Lmat2 matrix assembly equivalence test.
//
// Compares production, reference and direct particle assembly of the ECSIM
// response matrix on empty, dense, boundary and reassembly scenarios:
//   Reference: fill_matrixL2_Reference + stencil_Lmat2_Reference
//   Production: fill_matrixL2 + stencil_Lmat2
//   Direct particle assembly: direct per-particle rebuild from coordinates/physics
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
#include <stdexcept>
#include <string>
#include <vector>

#include "../unit/test_support.h"
#include "Mesh.h"
#include "Particle.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "bmatrix.h"
#include "grid/domain.h"
#include "memory.h"
#include "nlohmann/json.hpp"
#include "row_block.h"
#include "types.h"

// =============================================================================
// configuration constants
// =============================================================================
constexpr int NCELL = 6;             // cells per dimension
constexpr double CS = 0.1;           // cell size
constexpr double DT = 0.1;           // time step
constexpr double DICT_TOL = 1e-16;   // block-entry tolerance (same as blockToTriplets TOL)
constexpr double ZERO_TOL = 1e-14;   // structure filter: |v| <= ZERO_TOL is "zero"
constexpr double REL_TOL = 1e-9;     // relative value tolerance
constexpr double ABS_TOL = 1e-14;    // absolute value tolerance
constexpr int MAX_MISMATCH_PRINT = 10;
constexpr double SANITY_MIN_ABS_VALUE = 1e-6;   // sanity: matrix values must be non-trivial

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
            if (!std::isfinite(it.value()))
                throw std::runtime_error("nonfinite matrix entry before structure filtering");
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
// the same trusted accumulation path the reference assembly uses. Duplicate
// (row, col) triplets are summed by Eigen, matching the dict accumulation.
static Operator dictToCSR(const Dict& dict, int n) {
    std::vector<Trip> trips;
    trips.reserve(dictSize(dict));
    for (const auto& [row, cols] : dict) {
        for (const auto& [col, val] : cols) {
            if (!std::isfinite(val))
                throw std::runtime_error("nonfinite particle oracle entry");
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
// Direct particle response: independent of production block storage
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
static Dict assemble_response_directly_from_particles(const ParticlesArray& sp, const Field3d& fieldB,
                                                      const Domain& domain, double dt) {
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
            const double betaI = sp.mpw() * sp.charge / (1.0 + b.squaredNorm());
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
                                            if (!std::isfinite(val))
                                                throw std::runtime_error("nonfinite particle contribution");
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
    Operator reference_matrix;    // reference
    Operator production_matrix;   // production
    Dict particle_oracle;         // independent particle assembly
    long long reference_nnz = 0;
};

static nlohmann::json makeSpeciesConfig() {
    nlohmann::json cfg;
    cfg["Name"] = "lmat2_test_species";
    cfg["Charge"] = -1.0;
    cfg["Density"] = 1.0;
    cfg["Mass"] = 1.0;
    cfg["NumPartPerCell"] = 10;
    cfg["distribution"] = nlohmann::json::array();
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

// Runs all 3 assembly paths on a fresh/clean mesh state and compares them.
// Requires LmatX2 already prepared (by caller) and rowBlocksGlobal resized(0).
static void compare_production_reference_and_particle_response(
    const std::string& scenario, const std::string& stage, Mesh& mesh, ParticlesArray& sp, const Field3d& fieldB,
    const Domain& domain, SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>& rowBlocksGlobal, AssemblyResult& out) {
    const int n = domain.total_size() * 3;

    // reference: reference assembly. Fills the mesh.LmatX2 blocks, then assembles them.
    sp.fill_matrixL2_Reference<ShapeType::Linear>(mesh, fieldB, domain, DT);
    out.reference_matrix.resize(n, n);
    mesh.stencil_Lmat2_Reference(out.reference_matrix, domain);

    // production: current production assembly into per-thread row blocks, then merge.
    // rowBlocksGlobal must be reset before each use.
    for (int i = 0; i < static_cast<int>(rowBlocksGlobal.size()); ++i) {
        rowBlocksGlobal[i].resizeAndReset(0);
    }
    sp.fill_matrixL2<ShapeType::Linear>(mesh, fieldB, domain, DT, rowBlocksGlobal);
    out.production_matrix.resize(n, n);
    mesh.stencil_Lmat2(out.production_matrix, rowBlocksGlobal, mesh.workspacePtr);

    // independent particle assembly: direct per-particle dict oracle (independent of LmatX2 and
    // BlockDims indexing).
    out.particle_oracle = assemble_response_directly_from_particles(sp, fieldB, domain, DT);

    out.reference_nnz = extractEntries(out.reference_matrix, ZERO_TOL).size();

    std::cout << "  [" << stage << "] particles=" << sp.get_total_num_of_particles()
              << " nnz(reference)=" << out.reference_nnz << " nnz(particle oracle)=" << dictSize(out.particle_oracle)
              << std::endl;

    // -- sanity: inputs and results must be non-trivial ---------------------
    // The matrix value is betaL * s1*s2 * matB with
    //   betaL = 0.25*dt^2*(q/m)*(mpw*q/(1+b^2)); it depends on coord, B,
    //   q/m, mpw, dt (particle velocity does NOT enter Lmat2). If any input
    //   silently parsed to zero, all 3 paths would agree on an all-zero
    //   matrix and the comparisons would pass vacuously. Guard against that.
    check::require(sp.charge != 0.0, scenario + " [" + stage + "] sanity: charge != 0");
    check::require(sp.mass() != 0.0, scenario + " [" + stage + "] sanity: mass != 0");
    check::require(sp.mpw() > 0.0, scenario + " [" + stage + "] sanity: mpw > 0");
    double maxFieldB = 0.0;
    for (int i = 0; i < fieldB.sizes().x(); ++i)
        for (int j = 0; j < fieldB.sizes().y(); ++j)
            for (int k = 0; k < fieldB.sizes().z(); ++k)
                for (int d = 0; d < 3; ++d) maxFieldB = std::max(maxFieldB, std::fabs(fieldB(i, j, k, d)));
    check::require(maxFieldB > 0.0, scenario + " [" + stage + "] sanity: fieldB != 0");
    // non-empty scenarios must produce a non-empty matrix with non-trivial values
    if (sp.get_total_num_of_particles() > 0) {
        check::require(out.reference_nnz > 0, scenario + " [" + stage + "] sanity: nnz(reference) > 0");
        double maxAbs = 0.0;
        for (int row = 0; row < out.reference_matrix.outerSize(); ++row)
            for (Operator::InnerIterator it(out.reference_matrix, row); it; ++it)
                maxAbs = std::max(maxAbs, std::fabs(it.value()));
        check::require(maxAbs > SANITY_MIN_ABS_VALUE, scenario + " [" + stage + "] sanity: max|Lmat2| > " +
                                                          std::to_string(SANITY_MIN_ABS_VALUE) + " (got " +
                                                          std::to_string(maxAbs) + ")");
    }

    const std::string production_comparison = scenario + " [" + stage + "] reference == production";
    const std::string oracle_comparison = scenario + " [" + stage + "] reference == particle oracle";

    check::require(compareCSR(out.reference_matrix, out.production_matrix, "reference vs production"),
                   production_comparison);
    Operator particle_matrix = dictToCSR(out.particle_oracle, n);
    check::require(compareCSR(out.reference_matrix, particle_matrix, "reference vs direct particle assembly"),
                   oracle_comparison);
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

static void run_particle_configuration(const std::string& scenario, void (*populate)(ParticlesArray&, std::mt19937&)) {
    Domain domain = makeDomain();
    BoundaryConditionHandler bc_handler;   // empty handler: no BC modifications
    Mesh mesh;
    mesh.init(domain, DT, bc_handler);

    Field3d fieldB(domain.size(), 3);
    fieldB.setZero();
    fillUniformFieldB(fieldB);

    ParticlesArray sp(makeSpeciesConfig(), domain);
    std::mt19937 gen(20260601);
    populate(sp, gen);
    sp.update_cells(domain);   // groups particles into cells (required by both assembly paths)

    prepareLmatX2(mesh, sp);

    std::cout << "--- scenario: " << scenario << " ---" << std::endl;
    SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>> rowBlocksGlobal(omp_get_max_threads());

    AssemblyResult out;
    compare_production_reference_and_particle_response(scenario, "assemble-1", mesh, sp, fieldB, domain,
                                                       rowBlocksGlobal, out);

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
        // operators and rowBlocksGlobal buffers are reused as-is
        compare_production_reference_and_particle_response(scenario, "assemble-2", mesh, sp, fieldB, domain,
                                                           rowBlocksGlobal, out);
    }
}

// ---------------------------- populate functions ----------------------------
static void populate_empty_species(ParticlesArray& sp, std::mt19937&) {
    (void) sp;
}

static void populate_one_particle_at_cell_center(ParticlesArray& sp, std::mt19937&) {
    // one particle exactly at the center of an interior cell
    sp.add_particle(Particle(Vector3R(0.25, 0.25, 0.25), Vector3R(0.0, 0.0, 0.0)));
}

static void populate_one_particle_on_cell_edge(ParticlesArray& sp, std::mt19937&) {
    // one particle exactly on cell-boundary coordinates (x = 0.1 = one cell size)
    sp.add_particle(Particle(Vector3R(0.1, 0.2, 0.3), Vector3R(0.0, 0.0, 0.0)));
}

static void populate_dense_shape_subcell(ParticlesArray& sp, std::mt19937& gen) {
    // 95 particles in one shape subcell exercise the >=64 branch, batches
    // of 16 and all tail widths (8+4+2+1).
    for (int n = 0; n < 95; ++n) {
        const double dx = uniform(gen, -0.005, 0.005);
        const double dy = uniform(gen, -0.005, 0.005);
        const double dz = uniform(gen, -0.005, 0.005);
        sp.add_particle(Particle(Vector3R(0.23 + dx, 0.23 + dy, 0.23 + dz),
                                 Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
    }
}

static void populate_every_cell(ParticlesArray& sp, std::mt19937& gen) {
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
                    sp.add_particle(
                        Particle(Vector3R(cx + dx, cy + dy, cz + dz),
                                 Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
                }
            }
        }
    }
}

static void populate_random_cloud(ParticlesArray& sp, std::mt19937& gen) {
    // ~200 particles uniformly random over the domain [0, 0.6)^3
    const int count = 200;
    for (int n = 0; n < count; ++n) {
        sp.add_particle(Particle(Vector3R(uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6)),
                                 Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
    }
}

static void populate_domain_faces_and_corners(ParticlesArray& sp, std::mt19937& gen) {
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
                    sp.add_particle(Particle(
                        Vector3R(std::min(std::max(cx + dx, 0.0), 0.599), std::min(std::max(cy + dy, 0.0), 0.599),
                                 std::min(std::max(cz + dz, 0.0), 0.599)),
                        Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
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

static void populate_cloud_for_reassembly(ParticlesArray& sp, std::mt19937& gen) {
    // like random_scatter; the second assembly (after displacement) is handled
    // inside run_particle_configuration for the "reassembly" tag.
    for (int n = 0; n < 120; ++n) {
        sp.add_particle(Particle(Vector3R(uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6), uniform(gen, 0.0, 0.6)),
                                 Vector3R(uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1), uniform(gen, -0.1, 0.1))));
    }
}

// =============================================================================
// main
// =============================================================================
int main() {
    return check::run([] {
        std::cout << "========================================" << std::endl;
        std::cout << "Running Lmat2 Assembly Equivalence Tests" << std::endl;
        std::cout << "OMP_NUM_THREADS = " << omp_get_max_threads() << std::endl;
        std::cout << "========================================\n" << std::endl;

        run_particle_configuration("empty", populate_empty_species);
        run_particle_configuration("single_center", populate_one_particle_at_cell_center);
        run_particle_configuration("single_edge", populate_one_particle_on_cell_edge);
        run_particle_configuration("all_in_one_cell", populate_dense_shape_subcell);
        run_particle_configuration("full_fill", populate_every_cell);
        run_particle_configuration("random_scatter", populate_random_cloud);
        run_particle_configuration("boundaries", populate_domain_faces_and_corners);
        run_particle_configuration("reassembly", populate_cloud_for_reassembly);

        std::cout << "Lmat2 assembly checks passed\n";
    });
}
