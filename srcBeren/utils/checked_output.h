#pragma once

#include <fstream>
#include <stdexcept>
#include <string>

inline void flush_output(std::ofstream& file, const std::string& filename) {
    file.flush();
    if (!file)
        throw std::runtime_error("Cannot write " + filename);
}

inline void close_output(std::ofstream& file, const std::string& filename) {
    flush_output(file, filename);
    file.close();
    if (!file)
        throw std::runtime_error("Cannot close " + filename);
}
