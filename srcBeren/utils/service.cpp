// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "util.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

bool create_directory(const std::string& path) {
#ifdef _WIN32
    if (CreateDirectory(path.c_str(), NULL) == 0) {
        if (GetLastError() == ERROR_ALREADY_EXISTS) {
            std::cerr << "Directory " << path << " already exists"
                      << "\n";
            return false;
        } else {
            std::cerr << "Failed to create directory " << path << "\n";
            return false;
        }
    }
#else
    if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
        if (errno == EEXIST) {
            std::cerr << "Directory " << path << " already exists"
                      << "\n";
            return false;
        } else {
            std::cerr << "Failed to create directory " << path << "\n";
            return false;
        }
    }
#endif

    std::cerr << "Create directory " << path << " : SUCCESS"
              << "\n";
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
