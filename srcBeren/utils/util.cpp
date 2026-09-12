// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "util.h"

#include <omp.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

bool create_directory(const std::string& path) {
    const std::filesystem::path fsPath = path;
    if (std::filesystem::exists(fsPath)) {
        std::cerr << "Directory (or file) " << std::filesystem::absolute(fsPath) << " is already exists" << std::endl;
        return false;
    }
    const bool res = std::filesystem::create_directory(fsPath);
    if (!res) {
        std::cerr << "Failed to create directory " << std::filesystem::absolute(fsPath) << std::endl;
        return false;
    }

    std::cerr << "Create directory " << fsPath << " : SUCCESS" << "\n";
    return true;
}

std::vector<std::string> split_string(const std::string& s, const char delim) {
    std::vector<std::string> elems;
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

void raiseNumThreadException(int64_t obtainedNthr, int64_t desiredNthr, const std::source_location& location) {
    std::stringstream ss;
    ss << "Error: number of requested OMP threads does not coincide with number of desired threads; obtained: "
       << obtainedNthr << " but desired " << desiredNthr << "\n";
    ss << "In function: " << location.function_name() << "\n";
    ss << "Source location: " << location.file_name() << ":" << location.line() << ":" << location.column() << "\n";
    throw std::runtime_error(ss.str());
}
