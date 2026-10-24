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
    std::lock_guard<std::mutex> lock(log_mutex);
    log_file.open(path);
    log_file << "=== beren3d log ===\n\n";
}

inline void set_timestep(int ts) {
    std::lock_guard<std::mutex> lock(log_mutex);
    current_timestep = ts;
}

namespace detail {
inline void write(const std::string &prefix, const std::string &msg, bool flush) {
    std::lock_guard<std::mutex> lock(log_mutex);
    if (current_timestep > 0)
        log_file << "[" << current_timestep << "] " << prefix << msg << "\n";
    else
        log_file << prefix << msg << "\n";
    if (flush)
        log_file.flush();
    std::cout << prefix << msg << "\n";
    if (flush)
        std::cout.flush();
}
}   // namespace detail

inline void info(const std::string &msg) {
    detail::write("", msg, /*flush=*/false);
}
inline void warn(const std::string &msg) {
    detail::write("WARN: ", msg, /*flush=*/true);
}
inline void error(const std::string &msg) {
    detail::write("ERROR: ", msg, /*flush=*/true);
}

inline void close() {
    std::lock_guard<std::mutex> lock(log_mutex);
    log_file.flush();
    log_file.close();
}

}   // namespace logger

#endif
