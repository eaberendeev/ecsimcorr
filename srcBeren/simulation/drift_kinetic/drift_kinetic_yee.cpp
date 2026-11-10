#include "drift_kinetic_yee.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "simulation/eccapfim/algorithms_eccapfim.h"

namespace driftKinetic {
namespace {

double spline1(double s) {
    s = std::abs(s);
    return s <= 1.0 ? 1.0 - s : 0.0;
}

double spline2(double s) {
    s = std::abs(s);
    if (s <= 0.5)
        return 0.75 - s * s;
    if (s < 1.5)
        return 0.5 * (1.5 - s) * (1.5 - s);
    return 0.0;
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

    void setup(const TrackSegment& s) {
        double S2[3][2][2][kWidth];   // [axis][shift][end][node]
        for (int a = 0; a < 3; ++a) {
            const double node = s.cell[a];
            const double frac = (node - s.r0[a]) / s.dr[a];
            for (int k = 0; k < kWidth; ++k) {
                const double g = s.cell[a] - 2 + k;
                S2[a][0][0][k] = spline2(g - s.r0[a]);
                S2[a][0][1][k] = spline2(g - s.r1[a]);
                S2[a][1][0][k] = spline2(g + 0.5 - s.r0[a]);
                S2[a][1][1][k] = spline2(g + 0.5 - s.r1[a]);
                along_ii[a][k] = spline1(g + 0.5 - s.rh[a]);
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
            for (int m = 0; m < 4; ++m) {
                const auto& Sb = S2[b][m / 2];
                const auto& Sc = S2[c][m % 2];
                for (int kb = 0; kb < kWidth; ++kb)
                    for (int kc = 0; kc < kWidth; ++kc)
                        across[a][m][kb][kc] = (Sb[0][kb] * Sc[1][kc] + Sb[1][kb] * Sc[0][kc] +
                                                2 * (Sb[0][kb] * Sc[0][kc] + Sb[1][kb] * Sc[1][kc])) /
                                               6.0;
            }
        }
    }
};

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
// touched by an endpoint give t in {0, 1} and are merged by `unique`.
std::vector<TrackSegment> split_at_half_planes(const Vector3R& r0, const Vector3R& r1) {
    check_track(r0, r1);
    const Vector3R dr = r1 - r0;
    std::vector<double> t{0.0, 1.0};
    for (int a = 0; a < 3; ++a) {
        const int n0 = static_cast<int>(std::floor(r0[a] + 0.5));
        const int n1 = static_cast<int>(std::floor(r1[a] + 0.5));
        for (int n = std::min(n0, n1) + 1; n <= std::max(n0, n1); ++n)
            t.push_back(std::clamp((n - 0.5 - r0[a]) / dr[a], 0.0, 1.0));
    }
    std::sort(t.begin(), t.end());
    t.erase(std::unique(t.begin(), t.end()), t.end());

    std::vector<TrackSegment> segments;
    segments.reserve(t.size() - 1);
    Vector3R end = r0;
    for (size_t i = 1; i < t.size(); ++i) {
        TrackSegment s;
        s.r0 = end;
        s.r1 = (i + 1 < t.size()) ? Vector3R(r0 + t[i] * dr) : r1;
        end = s.r1;
        s.rh = 0.5 * (s.r0 + s.r1);
        s.dr = s.r1 - s.r0;
        s.length = s.dr.length();
        for (int a = 0; a < 3; ++a) s.cell[a] = static_cast<int>(std::floor(s.rh[a] + 0.5));
        segments.push_back(s);
    }
    return segments;
}

int YeeCoupling::owner(int physical_index, int axis) const {
    return algorithmsECCAPFIM::node_owner(physical_index, axis, domain_, bc_);
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
        for (int i = 0; i < 3; ++i) {
            const int ii = owner(first[0] + i, 0);
            for (int j = 0; j < 3; ++j) {
                const int jj = owner(first[1] + j, 1);
                for (int k = 0; k < 3; ++k) {
                    const int kk = owner(first[2] + k, 2);
                    result[c] += w[0][i] * w[1][j] * w[2][k] * level(magnetic, ii, jj, kk, theta)[c];
                }
            }
        }
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
        for (int i = 0; i < 3; ++i) {
            const int ii = owner(first[0] + i, 0);
            for (int j = 0; j < 3; ++j) {
                const int jj = owner(first[1] + j, 1);
                for (int k = 0; k < 3; ++k) {
                    const int kk = owner(first[2] + k, 2);
#pragma omp atomic update
                    magnetization(ii, jj, kk, c) += value * w[0][i] * w[1][j] * w[2][k];
                }
            }
        }
    }
}

Fields YeeCoupling::interpolate(const Field3d& electric, const MagneticLevels& magnetic, const Vector3R& start,
                                const Vector3R& end, StepInterval interval) const {
    const Vector3R a = domain_.to_cell_coordinates(start), b = domain_.to_cell_coordinates(end);
    const auto segments = split_at_half_planes(a, b);
    const double total = (b - a).length();

    // |B| and b are the half-sums of the orbit ends at their own time levels.
    const Vector3R B0 = magnetic_at(magnetic, start, interval.begin);
    const Vector3R B1 = magnetic_at(magnetic, end, interval.end);
    const double len0 = B0.length(), len1 = B1.length();
    if (!(len0 > 0.0) || !(len1 > 0.0) || !std::isfinite(len0) || !std::isfinite(len1))
        throw std::runtime_error("drift_kinetic: zero or non-finite magnetic field at a guiding centre");
    Fields fields;
    fields.lenB = 0.5 * (len0 + len1);
    Vector3R mean_b = 0.5 * (B0 / len0 + B1 / len1);
    // grad|B| = sum_c b_c grad B_c with b / |b|^2, since mean_b is not a unit vector.
    const Vector3R hat_b = mean_b / mean_b.squaredNorm();
    const double theta_mid = 0.5 * (interval.begin + interval.end);
    const Vector3R h = domain_.cell_size();

    Vector3R E(0, 0, 0), gradB(0, 0, 0), rotB(0, 0, 0);
    SplineTables t;
    for (const auto& s : segments) {
        t.setup(s);
        int idx[3][kWidth + 1];
        for (int d = 0; d < 3; ++d)
            for (int k = 0; k <= kWidth; ++k) idx[d][k] = owner(s.cell[d] - 2 + k, d);

        Vector3R Es(0, 0, 0), gradBs(0, 0, 0), rotBs(0, 0, 0);
        // Electric field: nodes cell-1 .. cell+1, the same weights as the current.
        for (int kx = 1; kx < kWidth; ++kx)
            for (int ky = 1; ky < kWidth; ++ky)
                for (int kz = 1; kz < kWidth; ++kz) {
                    const int i = idx[0][kx], j = idx[1][ky], k = idx[2][kz];
                    Es[0] += electric(i, j, k, 0) * t.along_ii[0][kx] * t.across[0][0][ky][kz];
                    Es[1] += electric(i, j, k, 1) * t.along_ii[1][ky] * t.across[1][0][kz][kx];
                    Es[2] += electric(i, j, k, 2) * t.along_ii[2][kz] * t.across[2][0][kx][ky];
                }
        // Forward differences dB[a] = (B(g + e_a) - B(g)) / h_a at the midpoint level.
        for (int kx = 0; kx < kWidth; ++kx)
            for (int ky = 0; ky < kWidth; ++ky)
                for (int kz = 0; kz < kWidth; ++kz) {
                    const int k3[3] = {kx, ky, kz};
                    const Vector3R B = level(magnetic, idx[0][kx], idx[1][ky], idx[2][kz], theta_mid);
                    const Vector3R dB[3] = {
                        (level(magnetic, idx[0][kx + 1], idx[1][ky], idx[2][kz], theta_mid) - B) / h[0],
                        (level(magnetic, idx[0][kx], idx[1][ky + 1], idx[2][kz], theta_mid) - B) / h[1],
                        (level(magnetic, idx[0][kx], idx[1][ky], idx[2][kz + 1], theta_mid) - B) / h[2]};
                    double W[3][3];
                    for (int p = 0; p < 3; ++p) {
                        const int q = (p + 1) % 3, r = (p + 2) % 3;
                        W[p][p] = t.along_ii[p][k3[p]] * t.across[p][3][k3[q]][k3[r]];
                        W[p][q] = t.along_ij[p][k3[p]] * t.across[p][1][k3[q]][k3[r]];
                        W[p][r] = t.along_ij[p][k3[p]] * t.across[p][2][k3[q]][k3[r]];
                    }
                    for (int p = 0; p < 3; ++p) {
                        const int q = (p + 1) % 3, r = (p + 2) % 3;
                        rotBs[p] += dB[q][r] * W[q][r] - dB[r][q] * W[r][q];
                        gradBs[p] += hat_b[0] * dB[p][0] * W[p][0] + hat_b[1] * dB[p][1] * W[p][1] +
                                     hat_b[2] * dB[p][2] * W[p][2];
                    }
                }
        const double dtau = total > 0 ? s.length / total : 1.0;
        E += dtau * Es;
        gradB += dtau * gradBs;
        rotB += dtau * rotBs;
    }
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
    const auto segments = split_at_half_planes(a, b);
    const double total = (b - a).length();

    const Vector3R B0 = magnetic_at(magnetic, start, interval.begin);
    const Vector3R B1 = magnetic_at(magnetic, end, interval.end);
    const Vector3R mean_b = 0.5 * (B0.normalized() + B1.normalized());
    const Vector3R hat_b = mean_b / mean_b.squaredNorm();
    deposit_magnetization(magnetization, start, hat_b, 0.5 * moment_weight);
    deposit_magnetization(magnetization, end, hat_b, 0.5 * moment_weight);

    SplineTables t;
    for (const auto& s : segments) {
        t.setup(s);
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
#pragma omp atomic update
                    current(i, j, k, 0) += jx;
#pragma omp atomic update
                    current(i, j, k, 1) += jy;
#pragma omp atomic update
                    current(i, j, k, 2) += jz;
                }
    }
}

}   // namespace driftKinetic
