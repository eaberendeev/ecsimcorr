#include "coupling/density.h"

#include "coupling/stencil.h"
#include "particles/particle_array.h"
#include "runtime/profiler.h"

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
    const Domain& domain = species.get_domain();
    const Grid& grid = domain.grid;
    const double weight = species.is_neutral() ? species.mpw() : species.mpw() * species.charge();
    constexpr int N = 2 * ShapeSize;
    assert(densityOnGrid.sizes() == grid.size() && densityOnGrid.nd() == 1);
    densityOnGrid.setZero();

#pragma omp parallel for schedule(dynamic, 64)
    for (int cell = 0; cell < species.size(); ++cell) {
        for (const Particle& particle : species(cell)) {
            const Vector3R coord = domain.to_cell_coordinates(particle.coord);
            const Stencil3D<N> stencil = make_deposition_stencil<ShapeFn, N>(coord);
            const Vector3I start = grid.to_storage_index(stencil.origin);
            for (int x = 0; x < N; ++x) {
                const double wx = stencil(x, Axis::X);
                for (int y = 0; y < N; ++y) {
                    const double wxy = wx * stencil(y, Axis::Y);
                    for (int z = 0; z < N; ++z) {
                        const double contribution = weight * wxy * stencil(z, Axis::Z);
                        if (contribution != 0.0) {
                            const Vector3I node = start + Vector3I(x, y, z);
                            assert(densityOnGrid.checkIndex(node.x(), node.y(), node.z(), 0));
#pragma omp atomic update
                            densityOnGrid(node.x(), node.y(), node.z(), 0) += contribution;
                        }
                    }
                }
            }
        }
    }
}

template <ShapeFunction ShapeFn, int ShapeSize>
void deposit_density_impl(const ParticlesArray& species, Field3d& densityOnGrid) {
    RECORD_TIMER;
    const Array3D<std::vector<Particle>>& cells = species.cells();
    const Domain& domain = species.get_domain();
    const Grid& grid = domain.grid;
    const double weight = species.is_neutral() ? species.mpw() : species.mpw() * species.charge();
    constexpr int N = 2 * ShapeSize;
    assert(cells.size() == grid.size());
    assert(densityOnGrid.sizes() == grid.size() && densityOnGrid.nd() == 1);
    densityOnGrid.setZero();

#pragma omp parallel
    {
        timer::flatTimer timer("OMP section");
        constexpr int blockWidth = 2;
        constexpr int ompGranularity = std::max(1, 64 / blockWidth / blockWidth / blockWidth);
        constexpr int bufferSize = N + blockWidth - 1;
        const Vector3I sizes = cells.size();

#pragma omp for collapse(3) schedule(dynamic, ompGranularity)
        for (int i = 0; i < sizes.x(); i += blockWidth) {
            for (int j = 0; j < sizes.y(); j += blockWidth) {
                for (int k = 0; k < sizes.z(); k += blockWidth) {
                    // Bin indices refer to storage. The shape envelope begins
                    // one logical node before the cell, regardless of halo width.
                    const Vector3I origin = deposition_origin(grid.to_node_index({i, j, k}));
                    Buffer3d<bufferSize> buffer;
                    buffer.zero();
                    bool empty = true;
                    for (int x = 0; x < std::min(blockWidth, sizes.x() - i); ++x)
                        for (int y = 0; y < std::min(blockWidth, sizes.y() - j); ++y)
                            for (int z = 0; z < std::min(blockWidth, sizes.z() - k); ++z) {
                                const Vector3I cell_origin = origin + Vector3I(x, y, z);
                                for (const Particle& particle : cells(i + x, j + y, k + z)) {
                                    empty = false;
                                    const Stencil3D<N> stencil = make_stencil_at<ShapeFn, N>(
                                        domain.to_cell_coordinates(particle.coord), cell_origin);
                                    for (int a = 0; a < N; ++a) {
                                        const double wx = stencil(a, Axis::X);
                                        for (int b = 0; b < N; ++b) {
                                            const double wxy = wx * stencil(b, Axis::Y);
                                            for (int c = 0; c < N; ++c)
                                                buffer(x + a, y + b, z + c, 0) += weight * wxy * stencil(c, Axis::Z);
                                        }
                                    }
                                }
                            }
                    if (empty)
                        continue;
                    const Vector3I start = grid.to_storage_index(origin);
                    const bool isBoundaryBlock =
                        start.x() < 0 || start.y() < 0 || start.z() < 0 || start.x() + bufferSize > grid.size(0) ||
                        start.y() + bufferSize > grid.size(1) || start.z() + bufferSize > grid.size(2);
                    for (int x = 0; x < bufferSize; ++x)
                        for (int y = 0; y < bufferSize; ++y)
                            for (int z = 0; z < bufferSize; ++z) {
                                const Vector3I node = start + Vector3I(x, y, z);
                                const double contribution = buffer(x, y, z, 0);
                                if (isBoundaryBlock) {
                                    if (contribution == 0.0)
                                        continue;
                                    assert(densityOnGrid.checkIndex(node.x(), node.y(), node.z(), 0));
                                }
#pragma omp atomic update
                                densityOnGrid(node.x(), node.y(), node.z(), 0) += contribution;
                            }
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
