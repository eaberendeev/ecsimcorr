// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef OUTPUT_UTIL_H_
#define OUTPUT_UTIL_H_

#include <string>
#include <vector>

#include "containers.h"

// format of field2D data storage.
// 2D arrays row-major for each dimension one after another
constexpr inline int output_field_index2d(int i, int j, int dim, int size1,
                                          int size2) {
    return dim * size1 * size2 + (j + size2 * i);
}

void output_field_plane(const Field3d& field, const int3& start,
                        const int3& end, int pos, int dim, int maxDim,
                        const std::string& filename,
                        const std::string& sNumber);
void output_array3d_plane(const Array3D<double>& field, const int3& sizes,
                          int pos, int dim, const std::string& filename,
                          const std::string& sNumber);
void output_array2d(const std::vector<float>& vector, int isize1, int isize2,
                    const std::string& filename);

void write_field_to_file(const std::string& dataName, const Field3d& field);
void read_field_from_file(const std::string& dataName, Field3d& field);

#endif
