#pragma once
#include <Shape.h>
#include <containers.h>

#include "config.h"
#include "timer.h"

// Плоский буфер для тока jx/jy/jz; индексируем по (n,m,k)
template <int SMAX>
struct CurrentBuffer {
    static_assert(SMAX > 0, "SMAX must be positive");
    alignas(64) std::array<double, 3 * SMAX * SMAX * SMAX> data_;

    static inline constexpr int idx(int n, int m, int k, int dim) noexcept {
        return (n * (SMAX * SMAX) + m * SMAX + k) * 3 + dim;
    }

    const double& operator()(int n, int m, int k, int dim) const {
        return data_[idx(n, m, k, dim)];
    }
    double& operator()(int n, int m, int k, int dim) {
        return data_[idx(n, m, k, dim)];
    }

    CurrentBuffer& operator+=(const CurrentBuffer& other) noexcept {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }
    inline void zero() noexcept {
        data_.fill(0.0);
    }
};

// Перенос локального буфера в глобальное поле fieldJ (атомарно)
template <int SMAX>
inline void flush_current_buffer(Field3d& fieldJ, const CurrentBuffer<SMAX>& buf, int start_x, int start_y,
                                 int start_z) noexcept {
    RECORD_TIMER;

    for (int n = 0; n < SMAX; ++n) {
        const int ix = start_x + n + GHOST_CELLS;
        for (int m = 0; m < SMAX; ++m) {
            const int iy = start_y + m + GHOST_CELLS;
            for (int k = 0; k < SMAX; ++k) {
                const int iz = start_z + k + GHOST_CELLS;
#pragma omp atomic update
                fieldJ(ix, iy, iz, Axis::X) += buf(n, m, k, Axis::X);
#pragma omp atomic update
                fieldJ(ix, iy, iz, Axis::Y) += buf(n, m, k, Axis::Y);
#pragma omp atomic update
                fieldJ(ix, iy, iz, Axis::Z) += buf(n, m, k, Axis::Z);
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void decompose_current(const ParticleShape<ShapeFn, ShapeSize>& no, const ParticleShape<ShapeFn, ShapeSize>& sh,
                       const Vector3R& value, CurrentBuffer<ShapeSize + 1>& curBuf) {
    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    // from -1 to ShapeSize
    for (int n = 0; n < ShapeSize; ++n) {
        int idx = n + 1;
        int idx05 = n + sh.start_.x() - no.start_.x() + 1;
        for (int m = 0; m < ShapeSize; ++m) {
            int idy = m + 1;
            int idy05 = m + sh.start_.y() - no.start_.y() + 1;
            for (int k = 0; k < ShapeSize; ++k) {
                int idz = k + 1;
                int idz05 = k + sh.start_.z() - no.start_.z() + 1;
                curBuf(idx05, idy, idz, X) += value.x() * sh(n, X) * no(m, Y) * no(k, Z);
                curBuf(idx, idy05, idz, Y) += value.y() * no(n, X) * sh(m, Y) * no(k, Z);
                curBuf(idx, idy, idz05, Z) += value.z() * no(n, X) * no(m, Y) * sh(k, Z);
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void decompose_esirkepov_current(const ParticleShape<ShapeFn, ShapeSize>& start,
                                 const ParticleShape<ShapeFn, ShapeSize>& end, const double qx, const double qy,
                                 const double qz, CurrentBuffer<ShapeSize>& curBuf) {
    RECORD_TIMER;

    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    curBuf.zero();

    for (int n = 0; n < ShapeSize; ++n) {
        double dsx = end(n, X) - start(n, X);
        for (int m = 0; m < ShapeSize; ++m) {
            double dsy = end(m, Y) - start(m, Y);
            for (int k = 0; k < ShapeSize; ++k) {
                double dsz = end(k, Z) - start(k, Z);

                // Предвычисляем веса для каждой компоненты
                double w_jx = (end(m, Y) * (2 * end(k, Z) + start(k, Z)) + start(m, Y) * (2 * start(k, Z) + end(k, Z)));
                double w_jy = (end(n, X) * (2 * end(k, Z) + start(k, Z)) + start(n, X) * (2 * start(k, Z) + end(k, Z)));
                double w_jz = (end(m, Y) * (2 * end(n, X) + start(n, X)) + start(m, Y) * (2 * start(n, X) + end(n, X)));

                if (n == 0)
                    curBuf(n, m, k, X) = -qx * dsx * w_jx;

                if (n > 0 && n < ShapeSize - 1) {
                    curBuf(n, m, k, X) = curBuf(n - 1, m, k, X) - qx * dsx * w_jx;
                }

                if (m == 0)
                    curBuf(n, m, k, Y) = -qy * dsy * w_jy;
                if (m > 0 && m < ShapeSize - 1) {
                    curBuf(n, m, k, Y) = curBuf(n, m - 1, k, Y) - qy * dsy * w_jy;
                }
                if (k == 0)
                    curBuf(n, m, k, Z) = -qz * dsz * w_jz;
                if (k > 0 && k < ShapeSize - 1) {
                    curBuf(n, m, k, Z) = curBuf(n, m, k - 1, Z) - qz * dsz * w_jz;
                }
            }
        }
    }
}
