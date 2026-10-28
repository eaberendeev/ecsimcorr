#pragma once

#include <string>
#include <vector>

#include "containers.h"

// format of field2D data storage.
// 2D arrays row-major for each dimension one after another
constexpr inline int output_field_index2d(int i, int j, int dim, int size1, int size2) {
    return dim * size1 * size2 + (j + size2 * i);
}

void output_field_plane(const Field3d& field, const Vector3I& start, const Vector3I& end, int pos, int dim, int maxDim,
                        const std::string& filename, const std::string& sNumber);
void output_array2d(const std::vector<float>& vector, int isize1, int isize2, const std::string& filename);
