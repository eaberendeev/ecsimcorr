#pragma once

#include <array>
#include <vector>

#include "World.h"
#include "containers.h"

using SliceIndices = std::array<std::vector<int>, 3>;
struct FieldSlice {
    int axis;
    int position;
    Field3d values;   // One node in the normal direction, including transverse ghosts.
};
using FieldSlices = std::vector<FieldSlice>;

SliceIndices slice_indices(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                           const Domain& domain);
void include_slice_indices(SliceIndices& target, const SliceIndices& source);
FieldSlices make_field_slices(const Vector3I& sizes, int components, const SliceIndices& positions);
void add_field_to_slices(FieldSlices& slices, const Field3d& field, const Field3d* extra = nullptr);
