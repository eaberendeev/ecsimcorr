#include "SliceFiles.h"

#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "indexing.h"
#include "util.h"

namespace {
std::pair<int, int> get_mesh_dimensions(const Vector3I& sizes, const int dim) {
    switch (dim) {
        case Axis::X:
            return {sizes.y(), sizes.z()};
        case Axis::Y:
            return {sizes.x(), sizes.z()};
        case Axis::Z:
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
    if (dim == Axis::X) {
        fullFilename += "_PlaneX_" + endFilename;
        for (int j = 0; j < size1; j++) {
            for (int k = 0; k < size2; k++) {
                for (int d = 0; d < maxDim; d++) {
                    int index = ind(pos + start.x(), j + start.y(), k + start.z(), d, extents.x(), extents.y(),
                                    extents.z(), maxDim);
                    vectorData[output_field_index2d(j, k, d, size1, size2)] = static_cast<float>(field.data()[index]);
                }
            }
        }
    } else if (dim == Axis::Y) {
        fullFilename += "_PlaneY_" + endFilename;
        for (int i = 0; i < size1; i++) {
            for (int k = 0; k < size2; k++) {
                for (int d = 0; d < maxDim; d++) {
                    int index = ind(i + start.x(), pos + start.y(), k + start.z(), d, extents.x(), extents.y(),
                                    extents.z(), maxDim);

                    vectorData[output_field_index2d(i, k, d, size1, size2)] = static_cast<float>(field.data()[index]);
                }
            }
        }
    } else if (dim == Axis::Z) {
        fullFilename += "_PlaneZ_" + endFilename;
        for (int i = 0; i < size1; i++) {
            for (int j = 0; j < size2; j++) {
                for (int d = 0; d < maxDim; d++) {
                    int index = ind(i + start.x(), j + start.y(), pos + start.z(), d, extents.x(), extents.y(),
                                    extents.z(), maxDim);
                    vectorData[output_field_index2d(i, j, d, size1, size2)] = static_cast<float>(field.data()[index]);
                }
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
    file.write((char*) &vector[0], vector.size() * sizeof(vector[0]));
}
