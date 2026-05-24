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

template <ShapeFunction ShapeFn, int ShapeSize>
void decompose_esirkepov_current_optimizedV9A(const ParticleShape<ShapeFn, ShapeSize>& start,
                                              const ParticleShape<ShapeFn, ShapeSize>& end, const double qx,
                                              const double qy, const double qz, CurrentBuffer<ShapeSize>& curBuf) {
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

                if (n == 0) {
                    curBuf(n, m, k, X) = -qx * dsx * w_jx;
                } else if (n < ShapeSize - 1) {
                    curBuf(n, m, k, X) = curBuf(n - 1, m, k, X) - qx * dsx * w_jx;
                }

                if (m == 0) {
                    curBuf(n, m, k, Y) = -qy * dsy * w_jy;
                } else if (m < ShapeSize - 1) {
                    curBuf(n, m, k, Y) = curBuf(n, m - 1, k, Y) - qy * dsy * w_jy;
                }
                if (k == 0) {
                    curBuf(n, m, k, Z) = -qz * dsz * w_jz;
                } else if (k < ShapeSize - 1) {
                    curBuf(n, m, k, Z) = curBuf(n, m, k - 1, Z) - qz * dsz * w_jz;
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void decompose_esirkepov_current_optimizedV9B(const ParticleShape<ShapeFn, ShapeSize>& start,
                                              const ParticleShape<ShapeFn, ShapeSize>& end, const double qx,
                                              const double qy, const double qz, CurrentBuffer<ShapeSize>& curBuf) {
    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    Eigen::Matrix<double, ShapeSize, ShapeSize> prevN;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevM;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevK;

    for (int n = 0; n < ShapeSize; ++n) {
        double dsx = end(n, X) - start(n, X);
        for (int m = 0; m < ShapeSize; ++m) {
            double dsy = end(m, Y) - start(m, Y);
            for (int k = 0; k < ShapeSize; ++k) {
                double dsz = end(k, Z) - start(k, Z);

                // Предвычисляем веса для каждой компоненты
                const double w_jx =
                    (end(m, Y) * (2 * end(k, Z) + start(k, Z)) + start(m, Y) * (2 * start(k, Z) + end(k, Z)));
                const double w_jy =
                    (end(n, X) * (2 * end(k, Z) + start(k, Z)) + start(n, X) * (2 * start(k, Z) + end(k, Z)));
                const double w_jz =
                    (end(m, Y) * (2 * end(n, X) + start(n, X)) + start(m, Y) * (2 * start(n, X) + end(n, X)));

                if (n == 0) {
                    prevN(m, k) = -qx * dsx * w_jx;
                    curBuf(n, m, k, X) += prevN(m, k);
                } else if (n < ShapeSize - 1) {
                    // const double tmp = prevN;
                    prevN(m, k) -= qx * dsx * w_jx;
                    curBuf(n, m, k, X) += prevN(m, k);
                }

                if (m == 0) {
                    prevM(n, k) = -qy * dsy * w_jy;
                    curBuf(n, m, k, Y) += prevM(n, k);
                } else if (m < ShapeSize - 1) {
                    prevM(n, k) -= qy * dsy * w_jy;
                    curBuf(n, m, k, Y) += prevM(n, k);
                }

                if (k == 0) {
                    prevK(m, n) = -qz * dsz * w_jz;
                    curBuf(n, m, k, Z) += prevK(m, n);
                } else if (k < ShapeSize - 1) {
                    prevK(m, n) -= qz * dsz * w_jz;
                    curBuf(n, m, k, Z) += prevK(m, n);
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void decompose_esirkepov_current_optimizedV9C(const ParticleShape<ShapeFn, ShapeSize>& start,
                                              const ParticleShape<ShapeFn, ShapeSize>& end, const double qx,
                                              const double qy, const double qz, CurrentBuffer<ShapeSize>& curBuf) {
    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    Eigen::Matrix<double, ShapeSize, ShapeSize> prevN;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevM;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevK;

    prevN.fill(0.0);
    prevM.fill(0.0);
    prevK.fill(0.0);

    for (int n = 0; n < ShapeSize; ++n) {
        double dsx = end(n, X) - start(n, X);
        for (int m = 0; m < ShapeSize; ++m) {
            double dsy = end(m, Y) - start(m, Y);
            for (int k = 0; k < ShapeSize; ++k) {
                double dsz = end(k, Z) - start(k, Z);

                // Предвычисляем веса для каждой компоненты
                const double w_jx =
                    (end(m, Y) * (2 * end(k, Z) + start(k, Z)) + start(m, Y) * (2 * start(k, Z) + end(k, Z)));
                const double w_jy =
                    (end(n, X) * (2 * end(k, Z) + start(k, Z)) + start(n, X) * (2 * start(k, Z) + end(k, Z)));
                const double w_jz =
                    (end(m, Y) * (2 * end(n, X) + start(n, X)) + start(m, Y) * (2 * start(n, X) + end(n, X)));

                if (n < ShapeSize - 1) {
                    prevN(m, k) -= qx * dsx * w_jx;
                    curBuf(n, m, k, X) += prevN(m, k);
                }
                if (m < ShapeSize - 1) {
                    prevM(n, k) -= qy * dsy * w_jy;
                    curBuf(n, m, k, Y) += prevM(n, k);
                }
                if (k < ShapeSize - 1) {
                    prevK(m, n) -= qz * dsz * w_jz;
                    curBuf(n, m, k, Z) += prevK(m, n);
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void decompose_esirkepov_current_optimizedV10(const ParticleShape<ShapeFn, ShapeSize>& start,
                                              const ParticleShape<ShapeFn, ShapeSize>& end, const double qx,
                                              const double qy, const double qz, CurrentBuffer<ShapeSize>& curBuf) {
    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    Eigen::Matrix<double, ShapeSize, ShapeSize> prevN;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevM;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevK;

    prevN.fill(0.0);
    prevM.fill(0.0);
    prevK.fill(0.0);

    Eigen::Matrix<double, ShapeSize, ShapeSize> wx;
    Eigen::Matrix<double, ShapeSize, ShapeSize> wy;
    Eigen::Matrix<double, ShapeSize, ShapeSize> wz;

    Eigen::Vector<double, ShapeSize> dsx;
    Eigen::Vector<double, ShapeSize> dsy;
    Eigen::Vector<double, ShapeSize> dsz;

    // Предвычисляем веса для каждой компоненты
    for (int i = 0; i < ShapeSize; ++i) {
        for (int j = 0; j < ShapeSize; ++j) {
            wx(i, j) = end(i, Y) * (2 * end(j, Z) + start(j, Z)) + start(i, Y) * (2 * start(j, Z) + end(j, Z));
            wy(i, j) = end(i, X) * (2 * end(j, Z) + start(j, Z)) + start(i, X) * (2 * start(j, Z) + end(j, Z));
            wz(i, j) = end(i, Y) * (2 * end(j, X) + start(j, X)) + start(i, Y) * (2 * start(j, X) + end(j, X));
        }
    }
    for (int i = 0; i < ShapeSize; ++i) {
        dsx(i) = end(i, X) - start(i, X);
        dsy(i) = end(i, Y) - start(i, Y);
        dsz(i) = end(i, Z) - start(i, Z);
    }

    for (int n = 0; n < ShapeSize; ++n) {
        for (int m = 0; m < ShapeSize; ++m) {
            for (int k = 0; k < ShapeSize; ++k) {
                const double w_jx = wx(m, k);
                const double w_jy = wy(n, k);
                const double w_jz = wz(m, n);

                if (n < ShapeSize - 1) {
                    prevN(m, k) -= qx * dsx(n) * w_jx;
                    curBuf(n, m, k, X) += prevN(m, k);
                }
                if (m < ShapeSize - 1) {
                    prevM(n, k) -= qy * dsy(m) * w_jy;
                    curBuf(n, m, k, Y) += prevM(n, k);
                }
                if (k < ShapeSize - 1) {
                    prevK(m, n) -= qz * dsz(k) * w_jz;
                    curBuf(n, m, k, Z) += prevK(m, n);
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void decompose_esirkepov_current_optimizedV11(const ParticleShape<ShapeFn, ShapeSize>& start,
                                              const ParticleShape<ShapeFn, ShapeSize>& end, const double qx,
                                              const double qy, const double qz, CurrentBuffer<ShapeSize>& curBuf) {
    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    Eigen::Matrix<double, ShapeSize, ShapeSize> prevN;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevM;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevK;

    prevN.fill(0.0);
    prevM.fill(0.0);
    prevK.fill(0.0);

    Eigen::Matrix<double, ShapeSize, ShapeSize> wx;
    Eigen::Matrix<double, ShapeSize, ShapeSize> wy;
    Eigen::Matrix<double, ShapeSize, ShapeSize> wz;

    Eigen::Vector<double, ShapeSize> scaledDsx;
    Eigen::Vector<double, ShapeSize> scaledDsy;
    Eigen::Vector<double, ShapeSize> scaledDsz;

    // Предвычисляем веса для каждой компоненты
    for (int i = 0; i < ShapeSize; ++i) {
        for (int j = 0; j < ShapeSize; ++j) {
            wx(i, j) = end(i, Y) * (2.0 * end(j, Z) + start(j, Z)) + start(i, Y) * (2.0 * start(j, Z) + end(j, Z));
            wy(i, j) = end(i, X) * (2.0 * end(j, Z) + start(j, Z)) + start(i, X) * (2.0 * start(j, Z) + end(j, Z));
            wz(i, j) = end(i, Y) * (2.0 * end(j, X) + start(j, X)) + start(i, Y) * (2.0 * start(j, X) + end(j, X));
        }
    }
    for (int i = 0; i < ShapeSize; ++i) {
        scaledDsx(i) = qx * (end(i, X) - start(i, X));
        scaledDsy(i) = qy * (end(i, Y) - start(i, Y));
        scaledDsz(i) = qz * (end(i, Z) - start(i, Z));
    }

    for (int n = 0; n < ShapeSize; ++n) {
        for (int m = 0; m < ShapeSize; ++m) {
            for (int k = 0; k < ShapeSize; ++k) {
                if (n < ShapeSize - 1) {
                    prevN(m, k) -= scaledDsx(n) * wx(m, k);
                    curBuf(n, m, k, X) += prevN(m, k);
                }
                if (m < ShapeSize - 1) {
                    prevM(n, k) -= scaledDsy(m) * wy(n, k);
                    curBuf(n, m, k, Y) += prevM(n, k);
                }
                if (k < ShapeSize - 1) {
                    prevK(m, n) -= scaledDsz(k) * wz(m, n);
                    curBuf(n, m, k, Z) += prevK(m, n);
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize, int bufferSize>
void decompose_esirkepov_current_optimizedV12(const ParticleShape<ShapeFn, ShapeSize>& start,
                                              const ParticleShape<ShapeFn, ShapeSize>& end, const double qx,
                                              const double qy, const double qz, CurrentBuffer<bufferSize>& curBuf,
                                              int offsetX, int offsetY, int offsetZ) {
    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    Eigen::Matrix<double, ShapeSize, ShapeSize> prevN;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevM;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevK;

    prevN.fill(0.0);
    prevM.fill(0.0);
    prevK.fill(0.0);

    Eigen::Matrix<double, ShapeSize, ShapeSize> wx;
    Eigen::Matrix<double, ShapeSize, ShapeSize> wy;
    Eigen::Matrix<double, ShapeSize, ShapeSize> wz;

    Eigen::Vector<double, ShapeSize> scaledDsx;
    Eigen::Vector<double, ShapeSize> scaledDsy;
    Eigen::Vector<double, ShapeSize> scaledDsz;

    // Предвычисляем веса для каждой компоненты
    for (int i = 0; i < ShapeSize; ++i) {
        for (int j = 0; j < ShapeSize; ++j) {
            wx(i, j) = end(i, Y) * (2.0 * end(j, Z) + start(j, Z)) + start(i, Y) * (2.0 * start(j, Z) + end(j, Z));
            wy(i, j) = end(i, X) * (2.0 * end(j, Z) + start(j, Z)) + start(i, X) * (2.0 * start(j, Z) + end(j, Z));
            wz(i, j) = end(i, Y) * (2.0 * end(j, X) + start(j, X)) + start(i, Y) * (2.0 * start(j, X) + end(j, X));
        }
    }
    for (int i = 0; i < ShapeSize; ++i) {
        scaledDsx(i) = qx * (end(i, X) - start(i, X));
        scaledDsy(i) = qy * (end(i, Y) - start(i, Y));
        scaledDsz(i) = qz * (end(i, Z) - start(i, Z));
    }

    for (int n = 0; n < ShapeSize; ++n) {
        for (int m = 0; m < ShapeSize; ++m) {
            for (int k = 0; k < ShapeSize; ++k) {
                if (n < ShapeSize - 1) {
                    prevN(m, k) -= scaledDsx(n) * wx(m, k);
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, X) += prevN(m, k);
                }
                if (m < ShapeSize - 1) {
                    prevM(n, k) -= scaledDsy(m) * wy(n, k);
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, Y) += prevM(n, k);
                }
                if (k < ShapeSize - 1) {
                    prevK(m, n) -= scaledDsz(k) * wz(m, n);
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, Z) += prevK(m, n);
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize, int bufferSize>
void decompose_esirkepov_current(const ParticleShape<ShapeFn, ShapeSize>& start,
                                 const ParticleShape<ShapeFn, ShapeSize>& end, const double qx, const double qy,
                                 const double qz, CurrentBuffer<bufferSize>& curBuf, int offsetX, int offsetY,
                                 int offsetZ) {
    constexpr int X = Axis::X;
    constexpr int Y = Axis::Y;
    constexpr int Z = Axis::Z;

    Eigen::Matrix<double, ShapeSize, ShapeSize> prevN;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevM;
    Eigen::Matrix<double, ShapeSize, ShapeSize> prevK;

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

                if (n == 0) {
                    prevN(m, k) = -qx * dsx * w_jx;
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, X) += prevN(m, k);
                }
                if (n > 0 && n < ShapeSize - 1) {
                    // const double tmp = prevN;
                    prevN(m, k) -= qx * dsx * w_jx;
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, X) += prevN(m, k);
                }

                if (m == 0) {
                    prevM(n, k) = -qy * dsy * w_jy;
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, Y) += prevM(n, k);
                }
                if (m > 0 && m < ShapeSize - 1) {
                    prevM(n, k) -= qy * dsy * w_jy;
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, Y) += prevM(n, k);
                }

                if (k == 0) {
                    prevK(m, n) = -qz * dsz * w_jz;
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, Z) += prevK(m, n);
                }
                if (k > 0 && k < ShapeSize - 1) {
                    prevK(m, n) -= qz * dsz * w_jz;
                    curBuf(n + offsetX, m + offsetY, k + offsetZ, Z) += prevK(m, n);
                }
            }
        }
    }
}
