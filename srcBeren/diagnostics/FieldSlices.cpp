#include "FieldSlices.h"

#include <algorithm>

#include "util.h"

SliceIndices slice_indices(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                           const Domain& domain) {
    SliceIndices indices;
    const std::array<std::vector<double>, 3> coordinates{x, y, z};
    for (int axis = 0; axis < 3; ++axis)
        for (double coordinate : coordinates[axis])
            indices[axis].push_back(int_value(coordinate / domain.cell_size(axis)));
    return indices;
}

void include_slice_indices(SliceIndices& target, const SliceIndices& source) {
    for (int axis = 0; axis < 3; ++axis)
        for (int position : source[axis])
            if (std::find(target[axis].begin(), target[axis].end(), position) == target[axis].end())
                target[axis].push_back(position);
}

FieldSlices make_field_slices(const Vector3I& sizes, int components, const SliceIndices& positions) {
    FieldSlices slices;
    for (int axis = 0; axis < 3; ++axis)
        for (int position : positions[axis]) {
            if (position < 0 || position >= sizes[axis])
                continue;
            auto plane_size = sizes;
            plane_size[axis] = 1;
            slices.push_back({axis, position, Field3d(plane_size, components)});
            slices.back().values.setZero();
        }
    return slices;
}

void add_field_to_slices(FieldSlices& slices, const Field3d& field, const Field3d* extra) {
    for (auto& plane : slices) {
        const auto sizes = plane.values.sizes();
        for (int i = 0; i < sizes.x(); ++i)
            for (int j = 0; j < sizes.y(); ++j)
                for (int k = 0; k < sizes.z(); ++k) {
                    Vector3I node(i, j, k);
                    node[plane.axis] = plane.position;
                    for (int d = 0; d < field.nd(); ++d) {
                        double value = field(node.x(), node.y(), node.z(), d);
                        if (extra)
                            value += (*extra)(node.x(), node.y(), node.z(), d);
                        plane.values(i, j, k, d) += value;
                    }
                }
    }
}
