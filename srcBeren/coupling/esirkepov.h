#pragma once
#include <omp.h>

#include "coupling/stencil.h"
#include "grid/domain.h"
#include "include/config.h"
#include "particles/Particle.h"
#include "particles/particle_path.h"
#include "runtime/profiler.h"
#include "utils/containers.h"

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

    const double& operator()(int n, int m, int k, Axis axis) const {
        return (*this)(n, m, k, static_cast<int>(axis));
    }

    double& operator()(int n, int m, int k, Axis axis) {
        return (*this)(n, m, k, static_cast<int>(axis));
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

template <int SMAX>
inline void flush_current_buffer_with_check(Field3d& fieldJ, const CurrentBuffer<SMAX>& buf,
                                            const Vector3I& storage_origin) noexcept;

// Перенос локального буфера в глобальное поле fieldJ (атомарно)
template <int SMAX>
inline void flush_current_buffer(Field3d& fieldJ, const CurrentBuffer<SMAX>& buf,
                                 const Vector3I& storage_origin) noexcept {
    // A particle can start exactly on the upper wall. Its shape fits in ghost
    // nodes, but the enclosing buffer can include zero padding outside the grid.
    assert(fieldJ.nd() == 3);
    const Vector3I start = storage_origin;
    const Vector3I sizes = fieldJ.sizes();
    const bool fits = start.x() >= 0 && start.y() >= 0 && start.z() >= 0 && start.x() + SMAX <= sizes.x() &&
                      start.y() + SMAX <= sizes.y() && start.z() + SMAX <= sizes.z();
    if (!fits) {
        flush_current_buffer_with_check(fieldJ, buf, storage_origin);
        return;
    }
    for (int n = 0; n < SMAX; ++n) {
        const int ix = start.x() + n;
        for (int m = 0; m < SMAX; ++m) {
            const int iy = start.y() + m;
            for (int k = 0; k < SMAX; ++k) {
                const int iz = start.z() + k;
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
inline void flush_current_buffer_with_check(Field3d& fieldJ, const CurrentBuffer<SMAX>& buf,
                                            const Vector3I& storage_origin) noexcept {
    assert(fieldJ.nd() == 3);
    const Vector3I start = storage_origin;
    for (int n = 0; n < SMAX; ++n) {
        const int ix = start.x() + n;
        for (int m = 0; m < SMAX; ++m) {
            const int iy = start.y() + m;
            for (int k = 0; k < SMAX; ++k) {
                const int iz = start.z() + k;
                if (ix >= 0 && iy >= 0 && iz >= 0 && ix < fieldJ.sizes().x() && iy < fieldJ.sizes().y() &&
                    iz < fieldJ.sizes().z()) {
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

// Endpoint weights must share logical nodes so their difference is the local
// charge change. The buffer is rooted at that same origin.
template <int ShapeSize>
inline void decompose_esirkepov_current_reference(const Stencil3D<ShapeSize>& start, const Stencil3D<ShapeSize>& end,
                                                  const double qx, const double qy, const double qz,
                                                  CurrentBuffer<ShapeSize>& curBuf) {
    assert(start.origin == end.origin);
    constexpr Axis X = Axis::X;
    constexpr Axis Y = Axis::Y;
    constexpr Axis Z = Axis::Z;

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

template <int ShapeSize>
inline void decompose_esirkepov_current(const Stencil3D<ShapeSize>& start, const Stencil3D<ShapeSize>& end,
                                        const double qx, const double qy, const double qz,
                                        CurrentBuffer<ShapeSize>& curBuf) {
    assert(start.origin == end.origin);
    constexpr Axis X = Axis::X;
    constexpr Axis Y = Axis::Y;
    constexpr Axis Z = Axis::Z;

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

// Deposit one straight segment with |to[d] - from[d]| <= 1 in cell units.
// The timestep must enforce this displacement bound; neighboring cell indices
// alone do not imply it. Rebase the four-node stencil to contain both endpoint
// shapes. `to` is normalized; from_world is in world coordinates. The caller
// keeps the original dt normalization, including each leg of a reflected path.
template <ShapeFunction ShapeFn, int Size, typename Deposit>
void deposit_segment_with_local_stencil(const Vector3R& from_world, const Vector3R& to, const Domain& domain,
                                        const Deposit& deposit) {
    const Vector3R from = domain.to_cell_coordinates(from_world);
    const Vector3I origin = trajectory_stencil_origin<ShapeFn>(from, to);
    static_assert(Size == 4, "trajectory deposition uses a common four-node envelope");
    const Stencil3D<Size> first = make_stencil_at<ShapeFn, Size>(from, origin);
    const Stencil3D<Size> last = make_stencil_at<ShapeFn, Size>(to, origin);
    assert(stencil_support_fits<ShapeFn>(first, from));
    assert(stencil_support_fits<ShapeFn>(last, to));
    deposit(first, last, false);
}

// Keep the cell buffer for the common case. If an endpoint's shape leaves
// that stencil, rebase only that trajectory and scatter its local buffer.
// Both legs of a reflected path retain the full timestep normalization:
// charge differences cancel at the wall hit.
template <ShapeFunction ShapeFn, int Size, typename Deposit>
inline void visit_current_path(const Stencil3D<Size>& start_shape, const Vector3R& start, const Particle& trial,
                               const Domain& domain, const PathPolicy& policy, const Deposit& deposit) {
    static_assert(Size == 4, "trajectory deposition uses a common four-node envelope");
    ParticleCurrentPath path{trial.coord, trial.coord, false, {}};
    if (policy.active())
        path = current_particle_path(start, trial, domain, policy);
    if (policy.events && path.crossing.clipped)
        policy.events->record(trial, {path.crossing.fate, path.crossing.face});
    const Vector3R turn = domain.to_cell_coordinates(path.turn);
    const bool turn_fits = stencil_support_fits<ShapeFn>(start_shape, turn);
    Stencil3D<Size> turn_shape;
    if (turn_fits) {
        turn_shape = make_stencil_at<ShapeFn, Size>(turn, start_shape.origin);
        deposit(start_shape, turn_shape, true);
    } else {
        deposit_segment_with_local_stencil<ShapeFn, Size>(start, turn, domain, deposit);
    }
    if (path.reflected) {
        const Vector3R end = domain.to_cell_coordinates(path.end);
        if (turn_fits && stencil_support_fits<ShapeFn>(start_shape, end)) {
            const Stencil3D<Size> end_shape = make_stencil_at<ShapeFn, Size>(end, start_shape.origin);
            deposit(turn_shape, end_shape, true);
        } else {
            deposit_segment_with_local_stencil<ShapeFn, Size>(path.turn, end, domain, deposit);
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
inline void move_and_calc_current_impl_reference(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                                 double charge, double mpw, const double dt, Field3d& fieldJ,
                                                 const PathPolicy& policy = PathPolicy::none()) {
    RECORD_TIMER;

    constexpr int SMAX = 2 * ShapeSize;
    static_assert(SMAX == 4, "trajectory deposition uses a common four-node envelope");

    const double qx = charge * domain.cell_size().x() / (6 * dt) * mpw;
    const double qy = charge * domain.cell_size().y() / (6 * dt) * mpw;
    const double qz = charge * domain.cell_size().z() / (6 * dt) * mpw;

#pragma omp parallel for schedule(dynamic, 64)
    for (int pk = 0; pk < particlesData.capacity(); ++pk) {
        std::vector<Particle>& particles = particlesData(pk);
        if (particles.empty()) {
            continue;
        }

        const Vector3I cell_origin = deposition_origin(domain.to_cell_coordinates(particles.front().coord));
        Stencil3D<SMAX> start_shape;
        CurrentBuffer<SMAX> curBuf, cellBuf;
        cellBuf.zero();
        curBuf.zero();

        for (Particle& particle : particles) {
            Vector3R start = particle.coord;
            const Vector3R coord = domain.to_cell_coordinates(start);
            assert(deposition_origin(coord) == cell_origin);
            start_shape = make_deposition_stencil<ShapeFn, SMAX>(coord);
            particle.move(dt);

            visit_current_path<ShapeFn>(start_shape, start, particle, domain, policy,
                                        [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool cell_local) {
                                            decompose_esirkepov_current_reference(a, b, qx, qy, qz, curBuf);
                                            if (cell_local)
                                                cellBuf += curBuf;
                                            else
                                                flush_current_buffer(fieldJ, curBuf, to_storage_index(a.origin));
                                        });
        }
        flush_current_buffer(fieldJ, cellBuf, to_storage_index(cell_origin));
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
inline void move_and_calc_current_impl(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                       double charge, double mpw, const double dt, Field3d& fieldJ,
                                       const PathPolicy& policy = PathPolicy::none()) {
    RECORD_TIMER;

    constexpr int SMAX = 2 * ShapeSize;
    static_assert(SMAX == 4, "trajectory deposition uses a common four-node envelope");

    const double qx = charge * domain.cell_size().x() / (6 * dt) * mpw;
    const double qy = charge * domain.cell_size().y() / (6 * dt) * mpw;
    const double qz = charge * domain.cell_size().z() / (6 * dt) * mpw;

    const int sizeX = particlesData.size().x();
    const int sizeY = particlesData.size().y();
    const int sizeZ = particlesData.size().z();

#pragma omp parallel
    {
        timer::flatTimer timer("OMP section");

        // Outgoing particles can still occupy halo bins before boundary handling.
        // Visit every allocated bin, as the reference path does.
#pragma omp for schedule(dynamic, 64) collapse(3)
        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {
                    std::vector<Particle>& particles = particlesData(i, j, k);
                    if (particles.empty())
                        continue;

                    const Vector3I cell_origin = deposition_origin(domain.to_cell_coordinates(particles.front().coord));
                    CurrentBuffer<SMAX> cellBuf;
                    cellBuf.zero();

                    for (Particle& particle : particles) {
                        const Vector3R start = particle.coord;
                        assert(deposition_origin(domain.to_cell_coordinates(start)) == cell_origin);
                        const Stencil3D<SMAX> start_shape =
                            make_stencil_at<ShapeFn, SMAX>(domain.to_cell_coordinates(start), cell_origin);
                        particle.move(dt);

                        visit_current_path<ShapeFn>(
                            start_shape, start, particle, domain, policy,
                            [&](const Stencil3D<SMAX>& a, const Stencil3D<SMAX>& b, bool cell_local) {
                                if (cell_local) {
                                    decompose_esirkepov_current(a, b, qx, qy, qz, cellBuf);
                                } else {
                                    CurrentBuffer<SMAX> scratch;
                                    scratch.zero();
                                    decompose_esirkepov_current(a, b, qx, qy, qz, scratch);
                                    flush_current_buffer(fieldJ, scratch, to_storage_index(a.origin));
                                }
                            });
                    }
                    flush_current_buffer(fieldJ, cellBuf, to_storage_index(cell_origin));
                }
            }
        }
    }
}

template <ShapeType ST = SHAPE>
inline void move_and_calc_current(Array3D<std::vector<Particle>>& particlesData, const Domain& domain, double charge,
                                  double mpw, const double dt, Field3d& fieldJ,
                                  const PathPolicy& policy = PathPolicy::none()) {
    RECORD_TIMER;

    move_and_calc_current_impl<ShapeKernel<ST>::weight, ShapeKernel<ST>::deposit_half_width>(
        particlesData, domain, charge, mpw, dt, fieldJ, policy);
}

template <ShapeType ST = SHAPE>
inline void move_and_calc_current_reference(Array3D<std::vector<Particle>>& particlesData, const Domain& domain,
                                            double charge, double mpw, const double dt, Field3d& fieldJ,
                                            const PathPolicy& policy = PathPolicy::none()) {
    RECORD_TIMER;
    move_and_calc_current_impl_reference<ShapeKernel<ST>::weight, ShapeKernel<ST>::deposit_half_width>(
        particlesData, domain, charge, mpw, dt, fieldJ, policy);
}
