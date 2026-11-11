#include "yee_stencil_matrix.h"

#include <cmath>
#include <stdexcept>
#include <unsupported/Eigen/FFT>

void YeeStencilMatrix::init(const Eigen::Vector3i& first_node, const Eigen::Vector3i& nodes, const DofMap& dof,
                            const Vector& free_rows, const Eigen::SparseMatrix<double>& constant_part, bool periodic) {
    periodic_ = periodic;
    for (int a = 0; a < 3; ++a) nodes_[a] = nodes[a];
    for (int c1 = 0; c1 < 3; ++c1) {
        int base = 0;
        for (int c2 = 0; c2 < 3; ++c2) {
            Box& box = boxes_[c1][c2];
            box.base = base;
            for (int a = 0; a < 3; ++a) {
                const int lo = (a == c2 && a != c1) ? -2 : (a == c1 ? -1 : -2);
                const int hi = (a == c1 && a != c2) ? 2 : (a == c2 ? 1 : 2);
                box.lo[a] = lo;
                box.width[a] = hi - lo + 1;
            }
            base += box.width[0] * box.width[1] * box.width[2];
        }
        if (base != kSlots)
            throw std::logic_error("YeeStencilMatrix: unexpected slot count");
    }

    rows_ = 3 * nodes.prod();
    free_.assign(rows_, 0);
    for (int r = 0; r < rows_; ++r) free_[r] = free_rows[r] != 0;
    columns_.assign(static_cast<size_t>(rows_) * kSlots, 0);
    enabled_.assign(columns_.size(), 0);
#pragma omp parallel for schedule(static)
    for (int node = 0; node < nodes.prod(); ++node) {
        const int ix = node / (nodes.y() * nodes.z()), iy = (node / nodes.z()) % nodes.y(), iz = node % nodes.z();
        const Eigen::Vector3i storage = first_node + Eigen::Vector3i(ix, iy, iz);
        for (int c1 = 0; c1 < 3; ++c1) {
            const int row = 3 * node + c1;
            for (int c2 = 0; c2 < 3; ++c2) {
                const Box& box = boxes_[c1][c2];
                for (int ox = 0; ox < box.width[0]; ++ox)
                    for (int oy = 0; oy < box.width[1]; ++oy)
                        for (int oz = 0; oz < box.width[2]; ++oz) {
                            const size_t slot = static_cast<size_t>(row) * kSlots + box.base +
                                                (ox * box.width[1] + oy) * box.width[2] + oz;
                            const int col = dof(storage.x() + ox + box.lo[0], storage.y() + oy + box.lo[1],
                                                storage.z() + oz + box.lo[2], c2);
                            // Constrained rows keep only their identity diagonal.
                            const bool valid = col >= 0 && (free_[row] ? free_[col] != 0 : col == row);
                            columns_[slot] = valid ? col : row;
                            enabled_[slot] = valid;
                        }
            }
        }
    }

    disabled_slots_.clear();
    for (size_t slot = 0; slot < enabled_.size(); ++slot)
        if (!enabled_[slot])
            disabled_slots_.push_back(slot);

    constant_values_.assign(columns_.size(), 0.0);
    const Eigen::SparseMatrix<double, Eigen::RowMajor> rows_major = constant_part;
    for (int r = 0; r < rows_; ++r)
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(rows_major, r); it; ++it) {
            bool placed = false;
            for (int s = 0; s < kSlots && !placed; ++s) {
                const size_t slot = static_cast<size_t>(r) * kSlots + s;
                if (enabled_[slot] && columns_[slot] == it.col()) {
                    constant_values_[slot] += it.value();
                    placed = true;
                }
            }
            if (!placed && it.value() != 0.0)
                throw std::logic_error("YeeStencilMatrix: Maxwell operator outside the particle pattern");
        }

    block_slots_.assign(rows_, {});
    for (int r = 0; r < rows_; ++r)
        for (int s = 0; s < kSlots; ++s) {
            const size_t slot = static_cast<size_t>(r) * kSlots + s;
            if (enabled_[slot] && columns_[slot] / 3 == r / 3)
                block_slots_[r].emplace_back(static_cast<int>(slot), columns_[slot] % 3);
        }
    reset();
}

void YeeStencilMatrix::reset() {
    values_ = constant_values_;
}

void YeeStencilMatrix::multiply(const Vector& x, Vector& y) const {
    y.resize(rows_);
#pragma omp parallel for schedule(static)
    for (int r = 0; r < rows_; ++r) {
        const size_t first = static_cast<size_t>(r) * kSlots;
        double sum = 0.0;
        for (int s = 0; s < kSlots; ++s) sum += values_[first + s] * x[columns_[first + s]];
        y[r] = sum;
    }
}

void YeeStencilMatrix::factor(Preconditioner preconditioner) {
    if (preconditioner != Preconditioner::BlockJacobi && !periodic_)
        throw std::invalid_argument("YeeStencilMatrix: circulant preconditioner needs a periodic unconstrained grid");
    preconditioner_ = preconditioner;
    if (preconditioner != Preconditioner::Circulant)
        factor_block_jacobi();
    if (preconditioner != Preconditioner::BlockJacobi)
        factor_circulant();
}

void YeeStencilMatrix::apply_preconditioner(const Vector& x, Vector& y) const {
    switch (preconditioner_) {
        case Preconditioner::BlockJacobi:
            apply_block_jacobi(x, y);
            return;
        case Preconditioner::Circulant:
            apply_circulant(x, y);
            return;
        case Preconditioner::CirculantBlockJacobi: {
            // y = BJ x; y += C (x - A y); y += BJ (x - A y).
            Vector r(rows_), z(rows_), ay(rows_);
            apply_block_jacobi(x, y);
            multiply(y, ay);
            r = x - ay;
            apply_circulant(r, z);
            y += z;
            multiply(y, ay);
            r = x - ay;
            apply_block_jacobi(r, z);
            y += z;
            return;
        }
    }
}

// Symbol of the translation average: (A v)_r = sum_o K_o v_{r+o} has the
// Fourier symbol sum_o K_o exp(+i k.o) for v_k = sum_r v_r exp(-i k.r).
void YeeStencilMatrix::factor_circulant() {
    const int nx = nodes_[0], ny = nodes_[1], nz = nodes_[2], nodes = nx * ny * nz;
    // Mean value of every slot over the nodes, for each row component.
    std::vector<double> mean(3 * kSlots, 0.0);
    for (int r = 0; r < rows_; ++r)
        for (int s = 0; s < kSlots; ++s) mean[(r % 3) * kSlots + s] += values_[static_cast<size_t>(r) * kSlots + s];
    for (double& value : mean) value /= nodes;

    const auto phase = [](int n) {
        std::vector<std::complex<double>> e(n * 5);   // exp(2 pi i k o / n) for o = -2..2
        for (int k = 0; k < n; ++k)
            for (int o = -2; o <= 2; ++o) e[k * 5 + o + 2] = std::polar(1.0, 2.0 * M_PI * k * o / n);
        return e;
    };
    const auto ex = phase(nx), ey = phase(ny), ez = phase(nz);
    symbol_inverse_.resize(nodes);
#pragma omp parallel for schedule(static)
    for (int mode = 0; mode < nodes; ++mode) {
        const int kx = mode / (ny * nz), ky = (mode / nz) % ny, kz = mode % nz;
        Eigen::Matrix3cd symbol = Eigen::Matrix3cd::Zero();
        for (int c1 = 0; c1 < 3; ++c1)
            for (int c2 = 0; c2 < 3; ++c2) {
                const Box& box = boxes_[c1][c2];
                for (int ox = 0; ox < box.width[0]; ++ox)
                    for (int oy = 0; oy < box.width[1]; ++oy)
                        for (int oz = 0; oz < box.width[2]; ++oz) {
                            const double k =
                                mean[c1 * kSlots + box.base + (ox * box.width[1] + oy) * box.width[2] + oz];
                            if (k == 0.0)
                                continue;
                            symbol(c1, c2) += k * ex[kx * 5 + ox + box.lo[0] + 2] * ey[ky * 5 + oy + box.lo[1] + 2] *
                                              ez[kz * 5 + oz + box.lo[2] + 2];
                        }
            }
        symbol_inverse_[mode] = symbol.inverse();
    }
}

namespace {
// In-place 3D FFT of `data` (x slowest, z fastest) along each axis.
void fft3(std::vector<std::complex<double>>& data, int nx, int ny, int nz, bool inverse) {
    const int n[3] = {nx, ny, nz};
    const int stride[3] = {ny * nz, nz, 1};
    for (int axis = 0; axis < 3; ++axis) {
        const int lines = nx * ny * nz / n[axis];
#pragma omp parallel
        {
            Eigen::FFT<double> fft;
            std::vector<std::complex<double>> in(n[axis]), out(n[axis]);
#pragma omp for schedule(static)
            for (int line = 0; line < lines; ++line) {
                // Index of the line start with the transformed axis removed.
                int start = 0, rest = line;
                for (int a = 2; a >= 0; --a) {
                    if (a == axis)
                        continue;
                    start += (rest % n[a]) * stride[a];
                    rest /= n[a];
                }
                for (int i = 0; i < n[axis]; ++i) in[i] = data[start + i * stride[axis]];
                if (inverse)
                    fft.inv(out, in);
                else
                    fft.fwd(out, in);
                for (int i = 0; i < n[axis]; ++i) data[start + i * stride[axis]] = out[i];
            }
        }
    }
}
}   // namespace

void YeeStencilMatrix::apply_circulant(const Vector& x, Vector& y) const {
    const int nx = nodes_[0], ny = nodes_[1], nz = nodes_[2], nodes = nx * ny * nz;
    std::vector<std::complex<double>> v[3];
    for (int c = 0; c < 3; ++c) {
        v[c].resize(nodes);
        for (int node = 0; node < nodes; ++node) v[c][node] = x[3 * node + c];
        fft3(v[c], nx, ny, nz, false);
    }
#pragma omp parallel for schedule(static)
    for (int mode = 0; mode < nodes; ++mode) {
        const Eigen::Vector3cd u = symbol_inverse_[mode] * Eigen::Vector3cd(v[0][mode], v[1][mode], v[2][mode]);
        for (int c = 0; c < 3; ++c) v[c][mode] = u[c];
    }
    y.resize(rows_);
    for (int c = 0; c < 3; ++c) {
        fft3(v[c], nx, ny, nz, true);
        for (int node = 0; node < nodes; ++node) y[3 * node + c] = v[c][node].real();
    }
}

void YeeStencilMatrix::factor_block_jacobi() {
    const int nodes = rows_ / 3;
    block_inverse_.resize(nodes);
#pragma omp parallel for schedule(static)
    for (int node = 0; node < nodes; ++node) {
        Eigen::Matrix3d block = Eigen::Matrix3d::Zero();
        for (int c1 = 0; c1 < 3; ++c1)
            for (const auto& [slot, c2] : block_slots_[3 * node + c1]) block(c1, c2) += values_[slot];
        block_inverse_[node] =
            std::abs(block.determinant()) > 0 ? Eigen::Matrix3d(block.inverse()) : Eigen::Matrix3d::Identity();
    }
}

void YeeStencilMatrix::apply_block_jacobi(const Vector& x, Vector& y) const {
    y.resize(rows_);
#pragma omp parallel for schedule(static)
    for (int node = 0; node < rows_ / 3; ++node) y.segment<3>(3 * node) = block_inverse_[node] * x.segment<3>(3 * node);
}

// GMRES(m) on A P^{-1} u = b, x = P^{-1} u, P the node-block Jacobi matrix.
bool YeeStencilMatrix::solve_gmres(const Vector& rhs, Vector& x, double tolerance, int restart, int max_iterations,
                                   int& iterations) const {
    iterations = 0;
    x = Vector::Zero(rows_);
    const double target = tolerance * rhs.norm();
    if (rhs.norm() == 0.0)
        return true;
    Vector r = rhs, w(rows_), z(rows_);
    std::vector<Vector> basis(restart + 1, Vector(rows_));
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(restart + 1, restart);
    Eigen::VectorXd cs(restart), sn(restart), g(restart + 1);
    while (iterations < max_iterations) {
        double beta = r.norm();
        if (beta <= target)
            return true;
        basis[0] = r / beta;
        g.setZero();
        g[0] = beta;
        int k = 0;
        for (; k < restart && iterations < max_iterations; ++k, ++iterations) {
            apply_preconditioner(basis[k], z);
            multiply(z, w);
            for (int i = 0; i <= k; ++i) {
                h(i, k) = w.dot(basis[i]);
                w -= h(i, k) * basis[i];
            }
            h(k + 1, k) = w.norm();
            if (h(k + 1, k) > 0)
                basis[k + 1] = w / h(k + 1, k);
            for (int i = 0; i < k; ++i) {
                const double t = cs[i] * h(i, k) + sn[i] * h(i + 1, k);
                h(i + 1, k) = -sn[i] * h(i, k) + cs[i] * h(i + 1, k);
                h(i, k) = t;
            }
            const double d = std::hypot(h(k, k), h(k + 1, k));
            cs[k] = h(k, k) / d;
            sn[k] = h(k + 1, k) / d;
            h(k, k) = d;
            h(k + 1, k) = 0;
            g[k + 1] = -sn[k] * g[k];
            g[k] = cs[k] * g[k];
            if (std::abs(g[k + 1]) <= target) {
                ++k;
                ++iterations;
                break;
            }
        }
        const Eigen::VectorXd y = h.topLeftCorner(k, k).triangularView<Eigen::Upper>().solve(g.head(k));
        Vector u = Vector::Zero(rows_);
        for (int i = 0; i < k; ++i) u += y[i] * basis[i];
        apply_preconditioner(u, z);
        x += z;
        multiply(x, w);
        r = rhs - w;
        if (!x.allFinite())
            return false;
    }
    return r.norm() <= target;
}

Eigen::SparseMatrix<double> YeeStencilMatrix::to_sparse() const {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(values_.size() / 4);
    for (int r = 0; r < rows_; ++r)
        for (int s = 0; s < kSlots; ++s) {
            const size_t slot = static_cast<size_t>(r) * kSlots + s;
            if (enabled_[slot] && values_[slot] != 0.0)
                triplets.emplace_back(r, columns_[slot], values_[slot]);
        }
    Eigen::SparseMatrix<double> result(rows_, rows_);
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}
