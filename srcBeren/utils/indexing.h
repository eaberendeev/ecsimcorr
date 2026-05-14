#pragma once

#include "types.h"
#include "vector3.h"

inline constexpr int ind(int x, int y, int z, int c, [[maybe_unused]] int Nx, int Ny, int Nz, int Nc) {
    return (c + Nc * (z + Nz * (y + Ny * x)));
}

struct IndexRange {
    Vector3I start;
    Vector3I end;
};

template <int SIZE_X, int SIZE_Y, int SIZE_Z, int DIMS = 1>
struct Indexer {
    static constexpr int size_x = SIZE_X;
    static constexpr int size_y = SIZE_Y;
    static constexpr int size_z = SIZE_Z;
    static constexpr int dims = DIMS;

    static constexpr int calculate(int x, int y, int z) noexcept {
        return x * (SIZE_Y * SIZE_Z) + y * SIZE_Z + z;
    }

    static constexpr int calculate(int x, int y, int z, int d) noexcept {
        return d + DIMS * (x * SIZE_Y * SIZE_Z + y * SIZE_Z + z);
    }
};
