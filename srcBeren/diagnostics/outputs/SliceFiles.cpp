#include "SliceFiles.h"

#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/indexing.h"
#include "output_util.h"
#include "util.h"

namespace {
std::pair<int, int> get_mesh_dimensions(const Vector3I& sizes, const int dim) {
    switch (dim) {
        case static_cast<int>(Axis::X):
            return {sizes.y(), sizes.z()};
        case static_cast<int>(Axis::Y):
            return {sizes.x(), sizes.z()};
        case static_cast<int>(Axis::Z):
            return {sizes.x(), sizes.y()};
        default:
            throw std::invalid_argument("Wrong dimension " + std::to_string(dim));
    }
}

}   // namespace

// writing Field3D array to 2D slice on dimension 'dim'
// and index 'pos' in this dimension
void output_field_plane(const Field3d& field, const Vector3I& start, const Vector3I& end, int pos, int dim, int maxDim,
                        const std::string& filename, const std::string& sNumber) {
    // check that slicing index is correct
    if (pos < 0 || pos >= field.sizes()[dim]) {
        // to do cout
        return;
    }
    const auto extents = field.sizes();
    const Vector3I sizes = end - start;
    auto meshDims = get_mesh_dimensions(sizes, dim);
    int size1 = meshDims.first;
    int size2 = meshDims.second;

    std::string fullFilename = filename;
    std::string endFilename = to_string(pos, 3) + "_" + sNumber;
    std::vector<float> vectorData(3 * size1 * size2);
    fullFilename += std::string("_Plane") + "XYZ"[dim] + "_" + endFilename;
    for (int u = 0; u < size1; ++u) {
        for (int v = 0; v < size2; ++v) {
            const int i = dim == static_cast<int>(Axis::X) ? pos : u;
            const int j = dim == static_cast<int>(Axis::Y) ? pos : (dim == static_cast<int>(Axis::X) ? u : v);
            const int k = dim == static_cast<int>(Axis::Z) ? pos : v;
            for (int d = 0; d < maxDim; ++d) {
                const int index =
                    ind(i + start.x(), j + start.y(), k + start.z(), d, extents.x(), extents.y(), extents.z(), maxDim);
                vectorData[output_field_index2d(u, v, d, size1, size2)] = static_cast<float>(field.data()[index]);
            }
        }
    }
    output_array2d(vectorData, size1, size2, fullFilename);
}

void output_array2d(const std::vector<float>& vector, int isize1, int isize2, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    float size1 = static_cast<float>(isize1);
    float size2 = static_cast<float>(isize2);
    file.write((char*) &size1, sizeof(size1));
    file.write((char*) &size2, sizeof(size2));
    file.write(reinterpret_cast<const char*>(vector.data()), vector.size() * sizeof(float));
    close_output(file, filename);
}
