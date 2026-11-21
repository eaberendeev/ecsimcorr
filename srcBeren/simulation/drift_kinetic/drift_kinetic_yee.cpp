#include "drift_kinetic_yee.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "simulation/eccapfim/algorithms_eccapfim.h"

namespace driftKinetic {
namespace {

// Linear and quadratic B-splines. Both branches are evaluated and selected,
// so that loops over stencil nodes vectorize.
inline double spline1(double s) {
    s = std::abs(s);
    return s <= 1.0 ? 1.0 - s : 0.0;
}

inline double spline2(double s) {
    s = std::abs(s);
    const double inner = 0.75 - s * s;
    const double outer = 0.5 * (1.5 - s) * (1.5 - s);
    return s <= 0.5 ? inner : (s < 1.5 ? outer : 0.0);
}

constexpr int kWidth = 4;   // nodes cell-2 .. cell+1 of the gradient stencil

// Separable factors of the segment weights on the nodes cell-2 .. cell+1.
// Every weight is a product of one `along` and one `across` entry.
struct SplineTables {
    double along_ii[3][kWidth];   // S1 along axis a at g + 1/2: E_a, J_a and d_a B_a
    double along_ij[3][kWidth];   // S1 along axis a at g + 1, split at the node: d_a B_c, c != a
    // Esirkepov factor over the axes b = a + 1, c = a + 2 (cyclic), at the
    // half-step shifts m = 2 shift_b + shift_c: 0 for E and J, 1 for d_a B_b,
    // 2 for d_a B_c, 3 for d_a B_a.
    double across[3][4][kWidth][kWidth];

    // `derivatives = false` fills only the E and J weights (along_ii, across m = 0).
    void setup(const TrackSegment& s, bool derivatives = true) {
        double S2[3][2][2][kWidth];   // [axis][shift][end][node]
        for (int a = 0; a < 3; ++a) {
            const double node = s.cell[a];
            const double frac = (node - s.r0[a]) / s.dr[a];
            for (int k = 0; k < kWidth; ++k) {
                const double g = s.cell[a] - 2 + k;
                S2[a][0][0][k] = spline2(g - s.r0[a]);
                S2[a][0][1][k] = spline2(g - s.r1[a]);
                along_ii[a][k] = spline1(g + 0.5 - s.rh[a]);
                if (!derivatives)
                    continue;
                S2[a][1][0][k] = spline2(g + 0.5 - s.r0[a]);
                S2[a][1][1][k] = spline2(g + 0.5 - s.r1[a]);
                // A derivative across a node plane inside the piece averages
                // the two half-pieces on either side of the node.
                if (0 < frac && frac < 1)
                    along_ij[a][k] = frac * spline1(g + 1.0 - 0.5 * (s.r0[a] + node)) +
                                     (1 - frac) * spline1(g + 1.0 - 0.5 * (s.r1[a] + node));
                else
                    along_ij[a][k] = spline1(g + 1.0 - s.rh[a]);
            }
        }
        for (int a = 0; a < 3; ++a) {
            const int b = (a + 1) % 3, c = (a + 2) % 3;
            for (int m = 0; m < (derivatives ? 4 : 1); ++m) {
                const auto& Sb = S2[b][m / 2];
                const auto& Sc = S2[c][m % 2];
                // (Sb0 (2 Sc0 + Sc1) + Sb1 (2 Sc1 + Sc0)) / 6: Esirkepov average over the piece.
                double u0[kWidth], u1[kWidth];
                for (int kc = 0; kc < kWidth; ++kc) {
                    u0[kc] = (2.0 * Sc[0][kc] + Sc[1][kc]) * (1.0 / 6.0);
                    u1[kc] = (2.0 * Sc[1][kc] + Sc[0][kc]) * (1.0 / 6.0);
                }
                for (int kb = 0; kb < kWidth; ++kb)
                    for (int kc = 0; kc < kWidth; ++kc) across[a][m][kb][kc] = Sb[0][kb] * u0[kc] + Sb[1][kb] * u1[kc];
            }
        }
    }
};

double component(const MagneticLevels& magnetic, int i, int j, int k, int c, double theta) {
    if (theta == 0.0)
        return (*magnetic.begin)(i, j, k, c);
    if (theta == 1.0)
        return (*magnetic.end)(i, j, k, c);
    if (theta == 0.5)
        return (*magnetic.half)(i, j, k, c);
    return (1 - theta) * (*magnetic.begin)(i, j, k, c) + theta * (*magnetic.end)(i, j, k, c);
}

Vector3R level(const MagneticLevels& magnetic, int i, int j, int k, double theta) {
    const auto at = [&](const Field3d& field) {
        return Vector3R(field(i, j, k, 0), field(i, j, k, 1), field(i, j, k, 2));
    };
    if (theta == 0.0)
        return at(*magnetic.begin);
    if (theta == 1.0)
        return at(*magnetic.end);
    if (theta == 0.5)
        return at(*magnetic.half);
    return (1 - theta) * at(*magnetic.begin) + theta * at(*magnetic.end);
}

void check_track(const Vector3R& a, const Vector3R& b) {
    for (int d = 0; d < 3; ++d)
        if (!std::isfinite(a[d]) || !std::isfinite(b[d]) || std::abs(a[d]) > 1e8 || std::abs(b[d]) > 1e8)
            throw std::runtime_error("drift_kinetic: invalid guiding-centre orbit coordinate");
}
}   // namespace

// Cell planes sit at n - 1/2; a piece takes the node of its midpoint. Planes
// touched by an endpoint give t in {0, 1} and are merged as duplicates.
// Short orbits use a stack buffer; a long orbit falls back to the heap.
template <typename F>
void for_each_segment(const Vector3R& r0, const Vector3R& r1, F&& f) {
    check_track(r0, r1);
    const Vector3R dr = r1 - r0;
    int n0[3], n1[3], count = 2;
    for (int a = 0; a < 3; ++a) {
        n0[a] = static_cast<int>(std::floor(r0[a] + 0.5));
        n1[a] = static_cast<int>(std::floor(r1[a] + 0.5));
        count += std::abs(n1[a] - n0[a]);
    }
    double stack[32];
    std::vector<double> heap;
    double* t = stack;
    if (count > 32) {
        heap.resize(count);
        t = heap.data();
    }
    int m = 0;
    t[m++] = 0.0;
    t[m++] = 1.0;
    for (int a = 0; a < 3; ++a)
        for (int n = std::min(n0[a], n1[a]) + 1; n <= std::max(n0[a], n1[a]); ++n)
            t[m++] = std::clamp((n - 0.5 - r0[a]) / dr[a], 0.0, 1.0);
    std::sort(t, t + m);
    m = static_cast<int>(std::unique(t, t + m) - t);
    Vector3R end = r0;
    for (int i = 1; i < m; ++i) {
        TrackSegment s;
        s.r0 = end;
        s.r1 = (i + 1 < m) ? Vector3R(r0 + t[i] * dr) : r1;
        end = s.r1;
        s.rh = 0.5 * (s.r0 + s.r1);
        s.dr = s.r1 - s.r0;
        s.length = s.dr.length();
        for (int a = 0; a < 3; ++a) s.cell[a] = static_cast<int>(std::floor(s.rh[a] + 0.5));
        f(s);
    }
}

std::vector<TrackSegment> split_at_half_planes(const Vector3R& r0, const Vector3R& r1) {
    std::vector<TrackSegment> segments;
    for_each_segment(r0, r1, [&](const TrackSegment& s) { segments.push_back(s); });
    return segments;
}

YeeCoupling::YeeCoupling(const Domain& domain, const BoundaryConditionHandler* bc) : domain_(domain), bc_(bc) {
    for (int d = 0; d < 3; ++d) {
        periodic_[d] = !bc || bc->is_periodic(d);
        cells_[d] = domain.num_cells(d);
        storage_[d] = domain.size()[d];
    }
}

// Forward differences d_a B_c = (B_c(g + e_a) - B_c(g)) / h_a of the three
// levels at every owner node; component a * 3 + c. Valid until the levels change.
void YeeCoupling::prepare(const MagneticLevels& magnetic) {
    const Field3d* sources[3] = {magnetic.begin, magnetic.half, magnetic.end};
    const Vector3R h = domain_.cell_size();
    // Physical node range: the owners on periodic axes, every stored node
    // (ghosts past a wall included) on the others.
    int lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = periodic_[d] ? 0 : -GHOST_CELLS;
        hi[d] = periodic_[d] ? cells_[d] : storage_[d] - GHOST_CELLS;
    }
    for (int level = 0; level < 3; ++level) {
        const Field3d& B = *sources[level];
        Field3d& d = differences_[level];
        if (d.nd() != 9 || d.sizes() != domain_.size())
            d.resize(domain_.size(), 9);
#pragma omp parallel for schedule(static)
        for (int i = lo[0]; i < hi[0]; ++i)
            for (int j = lo[1]; j < hi[1]; ++j)
                for (int k = lo[2]; k < hi[2]; ++k) {
                    const int p[3] = {i, j, k};
                    int g[3], up[3];
                    bool has_up[3];
                    for (int a = 0; a < 3; ++a) {
                        g[a] = owner(p[a], a);
                        // The last stored node past a wall has no neighbour: zero difference.
                        has_up[a] = periodic_[a] || p[a] + 1 < hi[a];
                        up[a] = has_up[a] ? owner(p[a] + 1, a) : g[a];
                    }
                    for (int c = 0; c < 3; ++c) {
                        const double b = B(g[0], g[1], g[2], c);
                        d(g[0], g[1], g[2], c) = (B(up[0], g[1], g[2], c) - b) / h[0];
                        d(g[0], g[1], g[2], 3 + c) = (B(g[0], up[1], g[2], c) - b) / h[1];
                        d(g[0], g[1], g[2], 6 + c) = (B(g[0], g[1], up[2], c) - b) / h[2];
                    }
                }
    }
    prepared_ = magnetic;
}

Vector3R YeeCoupling::magnetic_at(const Field3d& magnetic, const Vector3R& position) const {
    const MagneticLevels levels{&magnetic, &magnetic, &magnetic};
    return magnetic_at(levels, position, 0.0);
}

Vector3R YeeCoupling::magnetic_at(const MagneticLevels& magnetic, const Vector3R& position, double theta) const {
    const Vector3R r = domain_.to_cell_coordinates(position);
    Vector3R result(0, 0, 0);
    for (int c = 0; c < 3; ++c) {
        int first[3];
        double w[3][3];
        for (int a = 0; a < 3; ++a) {
            const double x = r[a] - (a == c ? 0.0 : 0.5);
            first[a] = static_cast<int>(std::floor(x + 0.5)) - 1;
            for (int n = 0; n < 3; ++n) w[a][n] = spline2(x - (first[a] + n));
        }
        int idx[3][3];
        for (int a = 0; a < 3; ++a)
            for (int n = 0; n < 3; ++n) idx[a][n] = owner(first[a] + n, a);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    result[c] +=
                        w[0][i] * w[1][j] * w[2][k] * component(magnetic, idx[0][i], idx[1][j], idx[2][k], c, theta);
    }
    return result;
}

void YeeCoupling::deposit_magnetization(Field3d& magnetization, const Vector3R& position, const Vector3R& direction,
                                        double moment) const {
    const Vector3R r = domain_.to_cell_coordinates(position);
    for (int c = 0; c < 3; ++c) {
        int first[3];
        double w[3][3];
        for (int a = 0; a < 3; ++a) {
            const double x = r[a] - (a == c ? 0.0 : 0.5);
            first[a] = static_cast<int>(std::floor(x + 0.5)) - 1;
            for (int n = 0; n < 3; ++n) w[a][n] = spline2(x - (first[a] + n));
        }
        const double value = -moment * direction[c];
        int idx[3][3];
        for (int a = 0; a < 3; ++a)
            for (int n = 0; n < 3; ++n) idx[a][n] = owner(first[a] + n, a);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    magnetization(idx[0][i], idx[1][j], idx[2][k], c) += value * w[0][i] * w[1][j] * w[2][k];
    }
}

Fields YeeCoupling::interpolate(const Field3d& electric, const MagneticLevels& magnetic, const Vector3R& start,
                                const Vector3R& end, StepInterval interval, const Vector3R* start_magnetic) const {
    const bool prepared =
        prepared_.begin == magnetic.begin && prepared_.half == magnetic.half && prepared_.end == magnetic.end;
    return interpolate_impl(electric, magnetic, start, end, interval, start_magnetic, prepared);
}

Fields YeeCoupling::interpolate_reference(const Field3d& electric, const MagneticLevels& magnetic,
                                          const Vector3R& start, const Vector3R& end, StepInterval interval) const {
    return interpolate_impl(electric, magnetic, start, end, interval, nullptr, false);
}

Vector3R YeeCoupling::start_magnetic(const MagneticLevels& magnetic, const Vector3R& start,
                                     StepInterval interval) const {
    return magnetic_at(magnetic, start, interval.begin);
}

Fields YeeCoupling::interpolate_impl(const Field3d& electric, const MagneticLevels& magnetic, const Vector3R& start,
                                     const Vector3R& end, StepInterval interval, const Vector3R* start_magnetic,
                                     bool prepared) const {
    const Vector3R a = domain_.to_cell_coordinates(start), b = domain_.to_cell_coordinates(end);
    const double total = (b - a).length();

    // |B| and b are the half-sums of the orbit ends at their own time levels.
    const Vector3R B0 = start_magnetic ? *start_magnetic : magnetic_at(magnetic, start, interval.begin);
    const Vector3R B1 = magnetic_at(magnetic, end, interval.end);
    const double len0 = B0.length(), len1 = B1.length();
    if (!(len0 > 0.0) || !(len1 > 0.0) || !std::isfinite(len0) || !std::isfinite(len1))
        throw std::runtime_error("drift_kinetic: zero or non-finite magnetic field at a guiding centre");
    Fields fields;
    fields.lenB = 0.5 * (len0 + len1);
    Vector3R mean_b = 0.5 * (B0 / len0 + B1 / len1);
    // grad|B| = sum_c b_c grad B_c with b / |b|^2, since mean_b is not a unit vector.
    const Vector3R hat_b = mean_b / mean_b.squaredNorm();
    const double theta = 0.5 * (interval.begin + interval.end);
    const Vector3R h = domain_.cell_size();
    // Precomputed differences of the level B(theta); a blend of B^n and B^{n+1} otherwise.
    const Field3d* d_single = theta == 0.0 ? &differences_[0] : (theta == 1.0 ? &differences_[2] : &differences_[1]);
    const bool blend = theta != 0.0 && theta != 0.5 && theta != 1.0;

    Vector3R E(0, 0, 0), gradB(0, 0, 0), rotB(0, 0, 0);
    SplineTables t;
    for_each_segment(a, b, [&](const TrackSegment& s) {
        t.setup(s);
        int idx[3][kWidth + 1];
        for (int d = 0; d < 3; ++d)
            for (int k = 0; k <= kWidth; ++k) idx[d][k] = owner(s.cell[d] - 2 + k, d);

        Vector3R Es(0, 0, 0), gradBs(0, 0, 0), rotBs(0, 0, 0);
        // Electric field: nodes cell-1 .. cell+1, the same weights as the current.
        for (int kx = 1; kx < kWidth; ++kx)
            for (int ky = 1; ky < kWidth; ++ky)
                for (int kz = 1; kz < kWidth; ++kz) {
                    const double* e = &electric(idx[0][kx], idx[1][ky], idx[2][kz], 0);
                    Es[0] += e[0] * t.along_ii[0][kx] * t.across[0][0][ky][kz];
                    Es[1] += e[1] * t.along_ii[1][ky] * t.across[1][0][kz][kx];
                    Es[2] += e[2] * t.along_ii[2][kz] * t.across[2][0][kx][ky];
                }
        // S[a][c] = sum_g W_ac(g) d_a B_c(g): the segment average of each
        // derivative d_a B_c on nodes cell-2 .. cell+1, with its own Yee weight
        // (d_a B_a: m = 3; d_a B_{a+1}: m = 1; d_a B_{a+2}: m = 2).
        double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int kx = 0; kx < kWidth; ++kx)
            for (int ky = 0; ky < kWidth; ++ky)
                for (int kz = 0; kz < kWidth; ++kz) {
                    const int k3[3] = {kx, ky, kz};
                    const int i = idx[0][kx], j = idx[1][ky], k = idx[2][kz];
                    double dB[9];
                    if (prepared && !blend) {
                        const double* d = &(*d_single)(i, j, k, 0);
                        for (int q = 0; q < 9; ++q) dB[q] = d[q];
                    } else if (prepared) {
                        const double* d0 = &differences_[0](i, j, k, 0);
                        const double* d1 = &differences_[2](i, j, k, 0);
                        for (int q = 0; q < 9; ++q) dB[q] = (1 - theta) * d0[q] + theta * d1[q];
                    } else {
                        const Vector3R B = level(magnetic, i, j, k, theta);
                        const Vector3R up[3] = {level(magnetic, idx[0][kx + 1], j, k, theta),
                                                level(magnetic, i, idx[1][ky + 1], k, theta),
                                                level(magnetic, i, j, idx[2][kz + 1], theta)};
                        for (int p = 0; p < 3; ++p)
                            for (int c = 0; c < 3; ++c) dB[3 * p + c] = (up[p][c] - B[c]) / h[p];
                    }
                    for (int p = 0; p < 3; ++p) {
                        const int q = (p + 1) % 3, r = (p + 2) % 3;
                        const double along_ii = t.along_ii[p][k3[p]], along_ij = t.along_ij[p][k3[p]];
                        S[p][p] += dB[3 * p + p] * along_ii * t.across[p][3][k3[q]][k3[r]];
                        S[p][q] += dB[3 * p + q] * along_ij * t.across[p][1][k3[q]][k3[r]];
                        S[p][r] += dB[3 * p + r] * along_ij * t.across[p][2][k3[q]][k3[r]];
                    }
                }
        // curl B = (S[y][z] - S[z][y], ...), grad|B| = sum_c b_c S[a][c].
        for (int p = 0; p < 3; ++p) {
            const int q = (p + 1) % 3, r = (p + 2) % 3;
            rotBs[p] = S[q][r] - S[r][q];
            gradBs[p] = hat_b[0] * S[p][0] + hat_b[1] * S[p][1] + hat_b[2] * S[p][2];
        }
        const double dtau = total > 0 ? s.length / total : 1.0;
        E += dtau * Es;
        gradB += dtau * gradBs;
        rotB += dtau * rotBs;
    });
    fields.electric = E;
    fields.b = mean_b.normalized();
    fields.gradB = gradB;
    fields.rotB = rotB;
    return fields;
}

void YeeCoupling::deposit(Field3d& current, Field3d& magnetization, const MagneticLevels& magnetic,
                          const Vector3R& start, const Vector3R& end, const Vector3R& velocity, double charge_weight,
                          double moment_weight, StepInterval interval) const {
    const Vector3R a = domain_.to_cell_coordinates(start), b = domain_.to_cell_coordinates(end);
    const double total = (b - a).length();

    const Vector3R B0 = magnetic_at(magnetic, start, interval.begin);
    const Vector3R B1 = magnetic_at(magnetic, end, interval.end);
    const Vector3R mean_b = 0.5 * (B0.normalized() + B1.normalized());
    const Vector3R hat_b = mean_b / mean_b.squaredNorm();
    deposit_magnetization(magnetization, start, hat_b, 0.5 * moment_weight);
    deposit_magnetization(magnetization, end, hat_b, 0.5 * moment_weight);

    SplineTables t;
    for_each_segment(a, b, [&](const TrackSegment& s) {
        t.setup(s, false);
        int idx[3][kWidth];
        for (int d = 0; d < 3; ++d)
            for (int k = 0; k < kWidth; ++k) idx[d][k] = owner(s.cell[d] - 2 + k, d);
        const double q = charge_weight * (total > 0 ? s.length / total : 1.0);
        for (int kx = 1; kx < kWidth; ++kx)
            for (int ky = 1; ky < kWidth; ++ky)
                for (int kz = 1; kz < kWidth; ++kz) {
                    const int i = idx[0][kx], j = idx[1][ky], k = idx[2][kz];
                    const double jx = q * velocity[0] * t.along_ii[0][kx] * t.across[0][0][ky][kz];
                    const double jy = q * velocity[1] * t.along_ii[1][ky] * t.across[1][0][kz][kx];
                    const double jz = q * velocity[2] * t.along_ii[2][kz] * t.across[2][0][kx][ky];
                    current(i, j, k, 0) += jx;
                    current(i, j, k, 1) += jy;
                    current(i, j, k, 2) += jz;
                }
    });
}

}   // namespace driftKinetic
