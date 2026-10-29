#pragma once

#include <array>
#include <cstddef>

#include "core/Vector3.h"

inline constexpr int ind(int x, int y, int z, int c, [[maybe_unused]] int Nx, int Ny, int Nz, int Nc) {
    return (c + Nc * (z + Nz * (y + Ny * x)));
}

// Decode one coordinate in the same row-major layout used by ind().
template <std::size_t N>
inline int index_coordinate(int index, int axis, const std::array<int, N>& dimensions) {
    int stride = 1;
    for (std::size_t i = axis + 1; i < N; ++i) stride *= dimensions[i];
    return (index / stride) % dimensions[axis];
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
