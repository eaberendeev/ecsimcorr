// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "output_util.h"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "containers.h"

void flush_output(std::ofstream& file, const std::string& filename) {
    file.flush();
    if (!file)
        throw std::runtime_error("Cannot write " + filename);
}

void close_output(std::ofstream& file, const std::string& filename) {
    flush_output(file, filename);
    file.close();
    if (!file)
        throw std::runtime_error("Cannot close " + filename);
}

void write_field_to_file(const std::string& dataName, const Field3d& field) {
    std::ofstream file_bin(dataName, std::ios::out | std::ios::binary);
    const auto sizes = field.sizes();
    const int size_x = sizes.x();
    const int size_y = sizes.y();
    const int size_z = sizes.z();

    file_bin.write((char*) &size_x, sizeof(size_x));
    file_bin.write((char*) &size_y, sizeof(size_y));
    file_bin.write((char*) &size_z, sizeof(size_z));

    const int dataSize = field.capacity();
    std::vector<double> data3d(dataSize);
    for (auto i = 0; i < dataSize; i++) {
        data3d[i] = field(i);
    }

    file_bin.write((char*) &data3d[0], dataSize * sizeof(data3d[0]));
    close_output(file_bin, dataName);
}

void read_field_from_file(const std::string& dataName, Field3d& field) {
    std::ifstream file_bin(dataName, std::ios::in | std::ios::binary);
    if (!file_bin.is_open()) {
        throw std::runtime_error("Cannot open " + dataName + " for reading");
    }
    int size_x = 0, size_y = 0, size_z = 0;
    file_bin.read(reinterpret_cast<char*>(&size_x), sizeof(size_x));
    file_bin.read(reinterpret_cast<char*>(&size_y), sizeof(size_y));
    file_bin.read(reinterpret_cast<char*>(&size_z), sizeof(size_z));
    if (!file_bin)
        throw std::runtime_error("Truncated recovery field header: " + dataName);
    const Vector3I sizes = field.sizes();
    if (field.nd() <= 0 || size_x != sizes.x() || size_y != sizes.y() || size_z != sizes.z())
        throw std::runtime_error("Recovery field dimensions do not match: " + dataName);

    // Validate the complete payload before replacing any live field values.
    std::vector<double> data3d(field.capacity());
    file_bin.read(reinterpret_cast<char*>(data3d.data()), data3d.size() * sizeof(double));
    if (!file_bin)
        throw std::runtime_error("Truncated recovery field payload: " + dataName);
    for (double value : data3d)
        if (!std::isfinite(value))
            throw std::runtime_error("Non-finite recovery field value: " + dataName);
    for (size_t i = 0; i < data3d.size(); ++i) field(i) = data3d[i];
    file_bin.close();
}
