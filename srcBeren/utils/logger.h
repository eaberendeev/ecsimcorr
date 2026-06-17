#pragma once

#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>

namespace logger {

inline std::ofstream log_file;
inline std::mutex log_mutex;
inline int current_timestep = 0;

inline void init(const std::string &path = "beren3d.log") {
    log_file.open(path);
    if (log_file.is_open()) log_file << "=== beren3d log ===\n\n";
}

inline void set_timestep(int ts) { current_timestep = ts; }

inline void info(const std::string &msg) {
    std::lock_guard<std::mutex> lock(log_mutex);
    if (log_file.is_open()) {
        if (current_timestep > 0)
            log_file << "[" << current_timestep << "] " << msg << "\n";
        else
            log_file << msg << "\n";
        log_file.flush();
    }
}

inline void warn(const std::string &msg) {
    std::lock_guard<std::mutex> lock(log_mutex);
    auto &out = (log_file.is_open()) ? log_file : std::cerr;
    if (current_timestep > 0)
        out << "[" << current_timestep << "] WARN: " << msg << "\n";
    else
        out << "WARN: " << msg << "\n";
    if (log_file.is_open()) log_file.flush();
}

inline void error(const std::string &msg) {
    std::lock_guard<std::mutex> lock(log_mutex);
    if (log_file.is_open())
        log_file << "ERROR: " << msg << "\n" << std::flush;
    std::cerr << "ERROR: " << msg << "\n" << std::flush;
}

inline void close() {
    if (log_file.is_open()) log_file.close();
}

}   // namespace logger

#endif
