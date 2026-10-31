#include "coupling/density.h"

#include "ParticlesArray.h"
#include "coupling/shape.h"
#include "timer.h"

namespace coupling {
namespace {
template <int size>
struct Buffer3d {
    static constexpr int dims = 1;

    static_assert(size > 0, "SMAX must be positive");
    alignas(64) std::array<double, dims * size * size * size> data_;

    static inline constexpr int idx(int n, int m, int k, int dim) noexcept {
        assert(0 <= n && n < size);
        assert(0 <= m && m < size);
        assert(0 <= k && k < size);
        assert(0 <= dim && dim < dims);
        return (n * (size * size) + m * size + k) * dims + dim;
    }

    double operator()(int n, int m, int k, int dim) const {
        return data_[idx(n, m, k, dim)];
    }
    double& operator()(int n, int m, int k, int dim) {
        return data_[idx(n, m, k, dim)];
    }

    Buffer3d& operator+=(const Buffer3d& other) noexcept {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }
    inline void zero() noexcept {
        data_.fill(0.0);
    }
};

template <ShapeFunction ShapeFn, int ShapeSize>
void deposit_density_impl_reference(const ParticlesArray& species, Field3d& densityOnGrid) {
    RECORD_TIMER;

    const auto& particlesData = species.particlesData;
    const Domain& domain_ = species.get_domain();
    const double mpw_ = species.mpw();
    const double charge = species.charge;
    constexpr auto SMAX = 2 * ShapeSize;
    densityOnGrid.setZero();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto j = 0; j < species.size(); ++j) {
        alignas(64) double sx[SMAX];
        alignas(64) double sy[SMAX];
        alignas(64) double sz[SMAX];

        for (const auto& particle : particlesData(j)) {
            // Vectorizable coordinate calculations
            const double ix = particle.coord.x() / domain_.cell_size().x();
            const double iy = particle.coord.y() / domain_.cell_size().y();
            const double iz = particle.coord.z() / domain_.cell_size().z();

            const int xk = floor(ix);
            const int yk = floor(iy);
            const int zk = floor(iz);

// Vectorizable shape calculations
#pragma omp simd
            for (int n = 0; n < SMAX; ++n) {
                sx[n] = ShapeFn(-ix + double(xk - GHOST_CELLS + n));
                sy[n] = ShapeFn(-iy + double(yk - GHOST_CELLS + n));
                sz[n] = ShapeFn(-iz + double(zk - GHOST_CELLS + n));
            }

            const double weight = species.is_neutral() ? mpw_ : mpw_ * charge;

            // Density accumulation with loop unrolling
            // #pragma unroll
            for (int n = 0; n < SMAX; ++n) {
                const int indx = xk + n;
                const double sxw = sx[n];

                for (int m = 0; m < SMAX; ++m) {
                    const int indy = yk + m;
                    const double sxyw = sxw * sy[m];

                    for (int k = 0; k < SMAX; ++k) {
                        const int indz = zk + k;
#pragma omp atomic update
                        densityOnGrid(indx, indy, indz, 0) += weight * sxyw * sz[k];
                    }
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void deposit_density_impl(const ParticlesArray& species, Field3d& densityOnGrid) {
    RECORD_TIMER;

    const auto& particlesData = species.particlesData;
    const Domain& domain_ = species.get_domain();
    const double mpw_ = species.mpw();
    const double charge = species.charge;
    constexpr auto SMAX = 2 * ShapeSize;
    densityOnGrid.setZero();

#pragma omp parallel
    {
        timer::flatTimer timer("OMP section");

        // TODO: this parameter have to be choosen individually for each particles distribution case, less part. tp cell
        // -> bigger
        constexpr int innerStep = 2;
        constexpr int ompGranularity = std::max(1, 64 / innerStep / innerStep / innerStep);

        constexpr int bufferSize = SMAX + innerStep - 1;

        const int sizeX = particlesData.size().x();
        const int sizeY = particlesData.size().y();
        const int sizeZ = particlesData.size().z();

#pragma omp for collapse(3) schedule(dynamic, ompGranularity)
        for (int i1 = GHOST_CELLS; i1 < sizeX; i1 += innerStep) {
            for (int j1 = GHOST_CELLS; j1 < sizeY; j1 += innerStep) {
                for (int k1 = GHOST_CELLS; k1 < sizeZ; k1 += innerStep) {
                    Buffer3d<bufferSize> localBuffer;
                    localBuffer.zero();
                    bool isBufferEmpty = true;

                    for (int i2 = 0; i2 < std::min(innerStep, sizeX - i1); ++i2) {
                        for (int j2 = 0; j2 < std::min(innerStep, sizeY - j1); ++j2) {
                            for (int k2 = 0; k2 < std::min(innerStep, sizeZ - k1); ++k2) {
                                const int i = i1 + i2;
                                const int j = j1 + j2;
                                const int k = k1 + k2;

                                alignas(64) double sx[SMAX];
                                alignas(64) double sy[SMAX];
                                alignas(64) double sz[SMAX];

                                // std::cout<<"Start particles data block!"<<std::endl;

                                if (particlesData(i, j, k).size() != 0)
                                    isBufferEmpty = false;

                                for (const auto& particle : particlesData(i, j, k)) {
                                    // Vectorizable coordinate calculations
                                    const double ix = particle.coord.x() / domain_.cell_size().x();
                                    const double iy = particle.coord.y() / domain_.cell_size().y();
                                    const double iz = particle.coord.z() / domain_.cell_size().z();
// Vectorizable shape calculations
#pragma omp simd
                                    for (int n = 0; n < SMAX; ++n) {
                                        sx[n] = ShapeFn(-ix + double(i - 2 * GHOST_CELLS + n));
                                        sy[n] = ShapeFn(-iy + double(j - 2 * GHOST_CELLS + n));
                                        sz[n] = ShapeFn(-iz + double(k - 2 * GHOST_CELLS + n));
                                    }

                                    const double weight = species.is_neutral() ? mpw_ : mpw_ * charge;

                                    // Density accumulation with loop unrolling
                                    // #pragma unroll
                                    for (int n = 0; n < SMAX; ++n) {
                                        const double sxw = sx[n];
                                        for (int m = 0; m < SMAX; ++m) {
                                            const double sxyw = sxw * sy[m];
                                            for (int k = 0; k < SMAX; ++k) {
                                                localBuffer(i2 + n, j2 + m, k2 + k, 0) += weight * sxyw * sz[k];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (!isBufferEmpty) {
                        const bool isBoundary = innerStep > 1 && (i1 + ShapeSize + innerStep - GHOST_CELLS >= sizeX ||
                                                                  j1 + ShapeSize + innerStep - GHOST_CELLS >= sizeY ||
                                                                  k1 + ShapeSize + innerStep - GHOST_CELLS >= sizeZ);
                        for (int n = 0; n < bufferSize; ++n) {
                            const int indx = i1 - GHOST_CELLS + n;
                            for (int m = 0; m < bufferSize; ++m) {
                                const int indy = j1 - GHOST_CELLS + m;
                                for (int k = 0; k < bufferSize; ++k) {
                                    const int indz = k1 - GHOST_CELLS + k;
                                    if (isBoundary) {
                                        if (indx < densityOnGrid.sizes().x() && indy < densityOnGrid.sizes().y() &&
                                            indz < densityOnGrid.sizes().z()) {
#pragma omp atomic update
                                            densityOnGrid(indx, indy, indz, 0) += localBuffer(n, m, k, 0);
                                        }
                                    } else {
#pragma omp atomic update
                                        densityOnGrid(indx, indy, indz, 0) += localBuffer(n, m, k, 0);
                                    }
                                }
                            }
                        }
                    }

                    // std::cout<<"End particles data block!"<<std::endl;
                }
            }
        }
    }
}

}   // namespace

template <ShapeType ST>
void deposit_density(const ParticlesArray& species, Field3d& density) {
    deposit_density_impl<ShapeKernel<ST>::weight, ShapeKernel<ST>::deposit_half_width>(species, density);
}

template <ShapeType ST>
void deposit_density_reference(const ParticlesArray& species, Field3d& density) {
    deposit_density_impl_reference<ShapeKernel<ST>::weight, ShapeKernel<ST>::deposit_half_width>(species, density);
}

template void deposit_density<ShapeType::Linear>(const ParticlesArray&, Field3d&);
template void deposit_density<ShapeType::Quadratic>(const ParticlesArray&, Field3d&);
template void deposit_density_reference<ShapeType::Linear>(const ParticlesArray&, Field3d&);
template void deposit_density_reference<ShapeType::Quadratic>(const ParticlesArray&, Field3d&);
}   // namespace coupling
