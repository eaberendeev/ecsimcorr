#pragma once
#include <Shape.h>
#include <containers.h>
#include <omp.h>

#include "Particle.h"
#include "World.h"
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

// Перенос локального буфера в глобальное поле fieldJ (атомарно)
template <int SMAX>
inline void flush_current_buffer_with_check(Field3d& fieldJ, const CurrentBuffer<SMAX>& buf, int start_x, int start_y,
                                            int start_z) noexcept {
    RECORD_TIMER;

    for (int n = 0; n < SMAX; ++n) {
        const int ix = start_x + n + GHOST_CELLS;
        for (int m = 0; m < SMAX; ++m) {
            const int iy = start_y + m + GHOST_CELLS;
            for (int k = 0; k < SMAX; ++k) {
                const int iz = start_z + k + GHOST_CELLS;
                if (ix < fieldJ.sizes().x() && iy < fieldJ.sizes().y() && iz < fieldJ.sizes().z()) {
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
}
template <ShapeFunction ShapeFn, int ShapeSize>
inline void decompose_current(const ParticleShape<ShapeFn, ShapeSize>& no, const ParticleShape<ShapeFn, ShapeSize>& sh,
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
inline void decompose_esirkepov_current_reference(const ParticleShape<ShapeFn, ShapeSize>& start,
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
inline void decompose_esirkepov_current(const ParticleShape<ShapeFn, ShapeSize>& start,
                                        const ParticleShape<ShapeFn, ShapeSize>& end, const double qx, const double qy,
                                        const double qz, CurrentBuffer<ShapeSize>& curBuf) {
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

    Eigen::Vector<double, ShapeSize - 1> scaledDsx;
    Eigen::Vector<double, ShapeSize - 1> scaledDsy;
    Eigen::Vector<double, ShapeSize - 1> scaledDsz;

    // Предвычисляем веса для каждой компоненты
    for (int i = 0; i < ShapeSize; ++i) {
        for (int j = 0; j < ShapeSize; ++j) {
            wx(i, j) = end(i, Y) * (2.0 * end(j, Z) + start(j, Z)) + start(i, Y) * (2.0 * start(j, Z) + end(j, Z));
            wy(i, j) = end(i, X) * (2.0 * end(j, Z) + start(j, Z)) + start(i, X) * (2.0 * start(j, Z) + end(j, Z));
            wz(i, j) = end(i, Y) * (2.0 * end(j, X) + start(j, X)) + start(i, Y) * (2.0 * start(j, X) + end(j, X));
        }
    }

    for (int i = 0; i < ShapeSize - 1; ++i) {
        scaledDsx(i) = qx * (end(i, X) - start(i, X));
        scaledDsy(i) = qy * (end(i, Y) - start(i, Y));
        scaledDsz(i) = qz * (end(i, Z) - start(i, Z));
    }

    for (int n = 0; n < ShapeSize - 1; ++n) {
        // use optimized matrix-scalar operation
        prevN -= scaledDsx(n) * wx;
        for (int k = 0; k < ShapeSize; ++k) {
            for (int m = 0; m < ShapeSize; ++m) {
                curBuf(n, m, k, X) += prevN(m, k);
            }
        }
    }

    for (int m = 0; m < ShapeSize - 1; ++m) {
        prevM -= scaledDsy(m) * wy;
        for (int k = 0; k < ShapeSize; ++k) {
            for (int n = 0; n < ShapeSize; ++n) {
                curBuf(n, m, k, Y) += prevM(n, k);
            }
        }
    }

    for (int k = 0; k < ShapeSize - 1; ++k) {
        prevK -= scaledDsz(k) * wz;
        for (int n = 0; n < ShapeSize; ++n) {
            for (int m = 0; m < ShapeSize; ++m) {
                curBuf(n, m, k, Z) += prevK(m, n);
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize, int bufferSize>
inline void decompose_esirkepov_current_offsetted(const ParticleShape<ShapeFn, ShapeSize>& start,
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

    Eigen::Vector<double, ShapeSize - 1> scaledDsx;
    Eigen::Vector<double, ShapeSize - 1> scaledDsy;
    Eigen::Vector<double, ShapeSize - 1> scaledDsz;

    // Предвычисляем веса для каждой компоненты
    for (int i = 0; i < ShapeSize; ++i) {
        for (int j = 0; j < ShapeSize; ++j) {
            wx(i, j) = end(i, Y) * (2.0 * end(j, Z) + start(j, Z)) + start(i, Y) * (2.0 * start(j, Z) + end(j, Z));
            wy(i, j) = end(i, X) * (2.0 * end(j, Z) + start(j, Z)) + start(i, X) * (2.0 * start(j, Z) + end(j, Z));
            wz(i, j) = end(i, Y) * (2.0 * end(j, X) + start(j, X)) + start(i, Y) * (2.0 * start(j, X) + end(j, X));
        }
    }

    for (int i = 0; i < ShapeSize - 1; ++i) {
        scaledDsx(i) = qx * (end(i, X) - start(i, X));
        scaledDsy(i) = qy * (end(i, Y) - start(i, Y));
        scaledDsz(i) = qz * (end(i, Z) - start(i, Z));
    }

    for (int n = 0; n < ShapeSize - 1; ++n) {
        // use optimized matrix-scalar operation
        prevN -= scaledDsx(n) * wx;
        for (int k = 0; k < ShapeSize; ++k) {
            for (int m = 0; m < ShapeSize; ++m) {
                curBuf(n + offsetX, m + offsetY, k + offsetZ, X) += prevN(m, k);
            }
        }
    }

    for (int m = 0; m < ShapeSize - 1; ++m) {
        prevM -= scaledDsy(m) * wy;
        for (int k = 0; k < ShapeSize; ++k) {
            for (int n = 0; n < ShapeSize; ++n) {
                curBuf(n + offsetX, m + offsetY, k + offsetZ, Y) += prevM(n, k);
            }
        }
    }

    for (int k = 0; k < ShapeSize - 1; ++k) {
        prevK -= scaledDsz(k) * wz;
        for (int n = 0; n < ShapeSize; ++n) {
            for (int m = 0; m < ShapeSize; ++m) {
                curBuf(n + offsetX, m + offsetY, k + offsetZ, Z) += prevK(m, n);
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
inline void move_and_calc_current_impl_reference(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                                 double charge, double mpw, const double dt, Field3d& fieldJ) {
    RECORD_TIMER;

    constexpr auto SMAX = 2 * ShapeSize;

    const double qx = charge * domain.cell_size().x() / (6 * dt) * mpw;
    const double qy = charge * domain.cell_size().y() / (6 * dt) * mpw;
    const double qz = charge * domain.cell_size().z() / (6 * dt) * mpw;

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < particlesData.capacity(); ++pk) {
        auto& particles = particlesData(pk);
        if (particles.empty()) {
            continue;
        }

        ParticleShape<ShapeFn, SMAX> start_shape;
        ParticleShape<ShapeFn, SMAX> end_shape;
        CurrentBuffer<SMAX> curBuf, cellBuf;
        cellBuf.zero();
        curBuf.zero();
        start_shape.fill_zero();

        for (auto& particle : particles) {
            Vector3R start = particle.coord;
            start_shape.fill_from_normalized(domain.to_cell_coordinates(start), GHOST_CELLS);
            particle.move(dt);

            Vector3R end = particle.coord;
            end_shape.fill_from_normalized(domain.to_cell_coordinates(end), start_shape.base_, GHOST_CELLS);
            decompose_esirkepov_current_reference(start_shape, end_shape, qx, qy, qz, curBuf);

            cellBuf += curBuf;
        }
        auto [start_x, start_y, start_z] = start_shape.start_.split();
        flush_current_buffer(fieldJ, cellBuf, start_x, start_y, start_z);
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
inline void move_and_calc_current_impl(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                       double charge, double mpw, const double dt, Field3d& fieldJ) {
    RECORD_TIMER;

    constexpr auto SMAX = 2 * ShapeSize;

    const double qx = charge * domain.cell_size().x() / (6 * dt) * mpw;
    const double qy = charge * domain.cell_size().y() / (6 * dt) * mpw;
    const double qz = charge * domain.cell_size().z() / (6 * dt) * mpw;

    // TODO: this parameter have to be choosen individually for each particles distribution case, less part. tp cell
    // -> bigger
    constexpr int innerStep = 1;
    constexpr int bufferSize = SMAX + innerStep - 1;
    constexpr int ompGranularity = std::max(1, 64 / innerStep / innerStep / innerStep);

#pragma omp parallel
    {
        timer::flatTimer timer("OMP section");

        ParticleShape<ShapeFn, SMAX> start_shape;
        ParticleShape<ShapeFn, SMAX> end_shape;
        CurrentBuffer<bufferSize> cellBuf;

        const int sizeX = particlesData.size().x();
        const int sizeY = particlesData.size().y();
        const int sizeZ = particlesData.size().z();

        // if this statement does not hold, then for case i1/j1/k1 = GHOST_CELLS object `cellBuf` will refer to
        // elements outside res
        static_assert(GHOST_CELLS + 1 >= ShapeSize);

// we don't have any particles in ghost cells, so do not take them into account
#pragma omp for schedule(dynamic, ompGranularity) collapse(3)
        for (int i1 = GHOST_CELLS; i1 < sizeX; i1 += innerStep) {
            for (int j1 = GHOST_CELLS; j1 < sizeY; j1 += innerStep) {
                for (int k1 = GHOST_CELLS; k1 < sizeZ; k1 += innerStep) {
                    int startX = -1;
                    int startY = -1;
                    int startZ = -1;

                    bool isBufferEmpty = true;

                    cellBuf.zero();

                    for (int i2 = 0; i2 < std::min(innerStep, sizeX - i1); ++i2) {
                        for (int j2 = 0; j2 < std::min(innerStep, sizeY - j1); ++j2) {
                            for (int k2 = 0; k2 < std::min(innerStep, sizeZ - k1); ++k2) {
                                const int i = i1 + i2;
                                const int j = j1 + j2;
                                const int k = k1 + k2;

                                auto& particles = particlesData(i, j, k);
                                if (particles.empty()) {
                                    continue;
                                }

                                timer::flatTimer timerCell("timer cell", particles.size());

                                const Vector3R coord0 = domain.to_cell_coordinates(particles[0].coord);
                                // integer coordinates of cell
                                const Vector3I baseCoord{static_cast<int>(std::floor(coord0.x())),
                                                         static_cast<int>(std::floor(coord0.y())),
                                                         static_cast<int>(std::floor(coord0.z()))};

                                start_shape.fill_zero();

                                for (auto& particle : particles) {
                                    Vector3R start = particle.coord;
                                    start_shape.fill_from_normalized(domain.to_cell_coordinates(start), baseCoord,
                                                                     GHOST_CELLS);
                                    particle.move(dt);

                                    Vector3R end = particle.coord;
                                    end_shape.fill_from_normalized(domain.to_cell_coordinates(end), start_shape.base_,
                                                                   GHOST_CELLS);
                                    if constexpr (innerStep == 1) {
                                        decompose_esirkepov_current(start_shape, end_shape, qx, qy, qz, cellBuf);
                                    } else {
                                        decompose_esirkepov_current_offsetted(start_shape, end_shape, qx, qy, qz,
                                                                              cellBuf, i2, j2, k2);
                                    }
                                }
                                if (isBufferEmpty) {
                                    auto [start_x, start_y, start_z] = start_shape.start_.split();
                                    startX = start_x - i2;
                                    startY = start_y - j2;
                                    startZ = start_z - k2;
                                }
                                isBufferEmpty = false;
                            }
                        }
                    }

                    if (!isBufferEmpty) {
                        const bool isBoundary = innerStep > 1 && (i1 + ShapeSize + innerStep - GHOST_CELLS >= sizeX ||
                                                                  j1 + ShapeSize + innerStep - GHOST_CELLS >= sizeY ||
                                                                  k1 + ShapeSize + innerStep - GHOST_CELLS >= sizeZ);
                        if (isBoundary)
                            flush_current_buffer_with_check<bufferSize>(fieldJ, cellBuf, startX, startY, startZ);
                        else
                            flush_current_buffer<bufferSize>(fieldJ, cellBuf, startX, startY, startZ);
                    }
                }
            }
        }
    }
}

inline void move_and_calc_current(Array3D<std::vector<Particle>>& particlesData, const Domain& domain, double charge,
                                  double mpw, const double dt, Field3d& fieldJ, ShapeType type) {
    switch (type) {
        case ShapeType::NGP:
            std::cout << "Move and calc current for NGP is not supported\n" << std::endl;
            exit(-1);
        case ShapeType::Linear:
            move_and_calc_current_impl<Shape, 2>(particlesData, domain, charge, mpw, dt, fieldJ);
            break;
        case ShapeType::Quadratic:
            move_and_calc_current_impl<Shape2, 2>(particlesData, domain, charge, mpw, dt, fieldJ);
            break;
    }
}

inline void move_and_calc_current_reference(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                            double charge, double mpw, const double dt, Field3d& fieldJ,
                                            ShapeType type) {
    switch (type) {
        case ShapeType::NGP:
            std::cout << "Move and calc current for NGP is not supported\n" << std::endl;
            exit(-1);
        case ShapeType::Linear:
            move_and_calc_current_impl_reference<Shape, 2>(particlesData, domain, charge, mpw, dt, fieldJ);
            break;
        case ShapeType::Quadratic:
            move_and_calc_current_impl_reference<Shape2, 2>(particlesData, domain, charge, mpw, dt, fieldJ);
            break;
    }
}
