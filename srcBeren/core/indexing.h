#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "core/vector3.h"

inline constexpr int ind(int x, int y, int z, int c, [[maybe_unused]] int Nx, int Ny, int Nz, int Nc) {
    return (c + Nc * (z + Nz * (y + Ny * x)));
}

// Floor a cell coordinate and apply a logical offset before narrowing to int.
// Checking the combined result also protects stencil origins such as floor(u)-1.
inline int checked_floor_index(double coordinate, int offset = 0) {
    const double lower = std::floor(coordinate);
    [[maybe_unused]] const int64_t min = int64_t(std::numeric_limits<int>::min()) - offset;
    [[maybe_unused]] const int64_t max = int64_t(std::numeric_limits<int>::max()) - offset;
    assert(std::isfinite(lower) && lower >= min && lower <= max);
    return static_cast<int>(static_cast<int64_t>(lower) + offset);
}

inline Vector3I checked_floor_index(const Vector3R& coordinate, int offset = 0) {
    return {checked_floor_index(coordinate.x(), offset), checked_floor_index(coordinate.y(), offset),
            checked_floor_index(coordinate.z(), offset)};
}

// Translate integer indices using only the halo offset, independently of Grid.
inline Vector3I to_storage_index(const Vector3I& node, int ghost_cells) {
    assert(ghost_cells >= 0);
    for ([[maybe_unused]] Axis axis : {Axis::X, Axis::Y, Axis::Z})
        assert(node[axis] <= std::numeric_limits<int>::max() - ghost_cells);
    return node + Vector3I(ghost_cells);
}
inline Vector3I to_node_index(const Vector3I& storage, int ghost_cells) {
    assert(ghost_cells >= 0);
    for ([[maybe_unused]] Axis axis : {Axis::X, Axis::Y, Axis::Z})
        assert(storage[axis] >= std::numeric_limits<int>::min() + ghost_cells);
    return storage - Vector3I(ghost_cells);
}

// Decode one coordinate in the same row-major layout used by ind().
template <int64_t N>
inline int64_t index_coordinate(int64_t index, int64_t axis, const std::array<int64_t, N>& dimensions) {
    assert(index >= 0 && axis >= 0 && axis < N);
    int64_t stride = 1;
    for (int64_t i = axis + 1; i < N; ++i) {
        assert(dimensions[i] > 0 && stride <= std::numeric_limits<int64_t>::max() / dimensions[i]);
        stride *= dimensions[i];
    }
    assert(dimensions[axis] > 0);
    return (index / stride) % dimensions[axis];
}

template <int64_t N>
inline int64_t index_coordinate(int64_t index, Axis axis, const std::array<int64_t, N>& dimensions) {
    return index_coordinate<N>(index, static_cast<int>(axis), dimensions);
}

struct IndexRange {
    Vector3I start;
    Vector3I end;

    int size() const {
        int res = 1;
        for (int i = 0; i < 3; ++i) {
            assert(start[i] <= end[i]);
            res *= (end[i] - start[i]);
        }
        return res;
    }
};

template <int SIZE_X, int SIZE_Y, int SIZE_Z, int DIMS = 1>
struct Indexer {
    static constexpr int size_x = SIZE_X;
    static constexpr int size_y = SIZE_Y;
    static constexpr int size_z = SIZE_Z;
    static constexpr int dims = DIMS;

    static constexpr int calculate(int x, int y, int z) noexcept {
        assert(x >= 0 && y >= 0 && z >= 0);
        assert(x < size_x && y < size_y && z < size_z);
        return ind(x, y, z, 0, SIZE_X, SIZE_Y, SIZE_Z, 1);
    }

    static constexpr int calculate(int x, int y, int z, int d) noexcept {
        assert(x >= 0 && y >= 0 && z >= 0 && d >= 0);
        assert(x < size_x && y < size_y && z < size_z && d < dims);
        return ind(x, y, z, d, SIZE_X, SIZE_Y, SIZE_Z, DIMS);
    }
};
