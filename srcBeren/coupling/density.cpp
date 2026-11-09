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

    // Only zero stencil padding may be omitted when flushing a boundary block.
    bool padding_is_zero(const Vector3I& begin, const Vector3I& end) const {
        for (int x = 0; x < size; ++x)
            for (int y = 0; y < size; ++y)
                for (int z = 0; z < size; ++z)
                    if (x < begin.x() || x >= end.x() || y < begin.y() || y >= end.y() || z < begin.z() || z >= end.z())
                        if ((*this)(x, y, z, 0) != 0.0)
                            return false;
        return true;
    }

    inline void zero() noexcept {
        data_.fill(0.0);
    }
};

// Specialize the interior path to keep compile-time loop bounds and unrolling.
// Only boundary blocks need runtime clipping of zero padding.
template <bool Boundary, int Size>
void flush_density_buffer(Field3d& density, const Buffer3d<Size>& buffer, const Vector3I& start) {
    Vector3I begin(0), end(Size);
    if constexpr (Boundary) {
        const Vector3I sizes = density.sizes();
        begin = Vector3I(std::max(0, -start.x()), std::max(0, -start.y()), std::max(0, -start.z()));
        end = Vector3I(std::min(Size, sizes.x() - start.x()), std::min(Size, sizes.y() - start.y()),
                       std::min(Size, sizes.z() - start.z()));
        assert(buffer.padding_is_zero(begin, end));
    }
    for (int x = begin.x(); x < end.x(); ++x)
        for (int y = begin.y(); y < end.y(); ++y)
            for (int z = begin.z(); z < end.z(); ++z) {
                const Vector3I node = start + Vector3I(x, y, z);
#pragma omp atomic update
                density(node.x(), node.y(), node.z(), 0) += buffer(x, y, z, 0);
            }
}

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
            const Vector3I start = to_storage_index(stencil.origin);
            for (int x = 0; x < N; ++x) {
                const double wx = stencil(x, Axis::X);
                for (int y = 0; y < N; ++y) {
                    const double wxy = wx * stencil(y, Axis::Y);
                    for (int z = 0; z < N; ++z) {
                        const double contribution = weight * wxy * stencil(z, Axis::Z);
                        if (contribution != 0.0) {
                            const Vector3I node = start + Vector3I(x, y, z);
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
                    const Vector3I start = deposition_origin(Vector3I(i, j, k));
                    const Vector3I origin = to_node_index(start);
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
                    const bool isBoundaryBlock =
                        start.x() < 0 || start.y() < 0 || start.z() < 0 || start.x() + bufferSize > sizes.x() ||
                        start.y() + bufferSize > sizes.y() || start.z() + bufferSize > sizes.z();
                    if (isBoundaryBlock)
                        flush_density_buffer<true>(densityOnGrid, buffer, start);
                    else
                        flush_density_buffer<false>(densityOnGrid, buffer, start);
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
