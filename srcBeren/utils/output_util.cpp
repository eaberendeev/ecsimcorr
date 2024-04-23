// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "output_util.h"

#include <fstream>
#include <iostream>
#include <vector>

#include "containers.h"
#include "service.h"
#include "util.h"

std::pair<int, int> get_mesh_dimensions(const int3& sizes, const int dim) {
    switch (dim) {
        case Dim::X:
            return {sizes.y(), sizes.z()};
        case Dim::Y:
            return {sizes.x(), sizes.z()};
        case Dim::Z:
            return {sizes.x(), sizes.y()};
        default:
            std::cout << "Wrong dimension " << dim << std::endl;
            return {-1, -1};
    }
}

// writing Field3D array to 2D slice on dimension 'dim'
// and index 'pos' in this dimension
void output_field_plane(const Field3d& field, const int3& start,
                        const int3& end, int pos, int dim, int maxDim,
                        const std::string& filename,
                        const std::string& sNumber) {
    // check that slicing index is correct
    if (pos < 0 || pos > field.size()(dim)) {
        // to do cout
        return;
    }
    const auto extents = field.size();
    const int3 sizes = end - start;
    auto meshDims = get_mesh_dimensions(sizes, dim);
    int size1 = meshDims.first;
    int size2 = meshDims.second;

    std::string fullFilename = filename;
    std::string endFilename = to_string(pos, 3) + "_" + sNumber;
    std::vector<float> vectorData(3 * size1 * size2);
    if (dim == Dim::X) {
        fullFilename += "_PlaneX_" + endFilename;
        for (int j = 0; j < size1; j++) {
            for (int k = 0; k < size2; k++) {
                for (int d = 0; d < maxDim; d++) {
                    int index =
                        ind(pos + start.x(), j + start.y(), k + start.z(), d,
                            extents.x(), extents.y(), extents.z(), maxDim);
                    vectorData[output_field_index2d(j, k, d, size1, size2)] =
                        static_cast<float>(field.data()[index]);
                }
            }
        }
    } else if (dim == Dim::Y) {
        fullFilename += "_PlaneY_" + endFilename;
        for (int i = 0; i < size1; i++) {
            for (int k = 0; k < size2; k++) {
                for (int d = 0; d < maxDim; d++) {
                    int index =
                        ind(i + start.x(), pos + start.y(), k + start.z(), d,
                            extents.x(), extents.y(), extents.z(), maxDim);

                    vectorData[output_field_index2d(i, k, d, size1, size2)] =
                        static_cast<float>(field.data()[index]);
                }
            }
        }
    } else if (dim == Dim::Z) {
        fullFilename += "_PlaneZ_" + endFilename;
        for (int i = 0; i < size1; i++) {
            for (int j = 0; j < size2; j++) {
                for (int d = 0; d < maxDim; d++) {
                    int index =
                        ind(i + start.x(), j + start.y(), pos + start.z(), d,
                            extents.x(), extents.y(), extents.z(), maxDim);
                    vectorData[output_field_index2d(i, j, d, size1, size2)] =
                        static_cast<float>(field.data()[index]);
                }
            }
        }
    }
    output_array2d(vectorData, size1, size2, fullFilename);
}

// writing 3D double array to 2D slice on dimension 'dim' and
// index 'pos' in this dimension
void output_array3d_plane(const Array3D<double>& array3D, const int3& sizes,
                          int pos, int dim, const std::string& filename,
                          const std::string& sNumber) {
    // check that slicing index is correct
    if (pos < 0 || pos > array3D.size()(dim)) {
        return;
    }
    auto meshDims = get_mesh_dimensions(sizes, dim);
    int size1 = meshDims.first;
    int size2 = meshDims.second;
    std::string fullFilename = filename;
    std::string endFilename = to_string(pos, 3) + "_" + sNumber;

    Array2D<float> array2D(size1, size2);
    if (dim == Dim::X) {
        fullFilename += "_PlaneX_" + endFilename;
        for (int j = 0; j < size1; j++) {
            for (int k = 0; k < size2; k++) {
                array2D(j, k) = static_cast<float>(array3D(pos, j, k));
            }
        }
    } else if (dim == Dim::Y) {
        fullFilename += "_PlaneY_" + endFilename;
        for (int i = 0; i < size1; i++) {
            for (int k = 0; k < size2; k++) {
                array2D(i, k) = static_cast<float>(array3D(i, pos, k));
            }
        }
    } else if (dim == Dim::Z) {
        fullFilename += "_PlaneZ_" + endFilename;
        for (int i = 0; i < size1; i++) {
            for (int j = 0; j < size2; j++) {
                array2D(i, j) = static_cast<float>(array3D(i, j, pos));
            }
        }
    }
    output_array2d(array2D.data(), size1, size2, fullFilename);
}

void output_array2d(const std::vector<float>& vector, int isize1, int isize2,
                    const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    float size1 = static_cast<float>(isize1);
    float size2 = static_cast<float>(isize2);
    file.write((char*) &size1, sizeof(size1));
    file.write((char*) &size2, sizeof(size2));
    file.write((char*) &vector[0], vector.size() * sizeof(vector[0]));
}

void write_field_to_file(const std::string& dataName, const Field3d& field) {
    std::ofstream file_bin(dataName, std::ios::out | std::ios::binary);
    const int size_x = field.size().x();
    const int size_y = field.size().y();
    const int size_z = field.size().z();

    file_bin.write((char*) &size_x, sizeof(size_x));
    file_bin.write((char*) &size_y, sizeof(size_y));
    file_bin.write((char*) &size_z, sizeof(size_z));

    const int dataSize = field.capacity();
    std::vector<double> data3d(dataSize);
    for (auto i = 0; i < dataSize; i++) {
        data3d[i] = field(i);
    }

    file_bin.write((char*) &data3d[0], dataSize * sizeof(data3d[0]));
    file_bin.close();
}

void read_field_from_file(const std::string& dataName, Field3d& field) {
    std::ifstream file_bin(dataName, std::ios::in | std::ios::binary);
    int size_x, size_y, size_z;
    file_bin.read((char*) &size_x, sizeof(size_x));
    file_bin.read((char*) &size_y, sizeof(size_y));
    file_bin.read((char*) &size_z, sizeof(size_z));
    std::cout << size_x << " " << size_y << " " << size_z << "\n";
    int3 size = field.size();
    if (size_x != size.x() || size_y != size.y() || size_z != size.z()) {
        std::cout << "Invalid reading from recovery! Field size is invalid! \n";
        exit(0);
    }

    const int dataSize = field.capacity();
    std::vector<double> data3d(dataSize);
    file_bin.read((char*) &data3d[0], dataSize * sizeof(data3d[0]));
    for (auto i = 0; i < dataSize; i++) {
        field(i) = data3d[i];
    }
    file_bin.close();
}

