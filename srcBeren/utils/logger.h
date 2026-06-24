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
    log_file << "=== beren3d log ===\n\n";
}

inline void set_timestep(int ts) { current_timestep = ts; }

namespace detail {
inline void write(const std::string &prefix, const std::string &msg) {
    std::lock_guard<std::mutex> lock(log_mutex);
    if (current_timestep > 0)
        log_file << "[" << current_timestep << "] " << prefix << msg << "\n";
    else
        log_file << prefix << msg << "\n";
    log_file.flush();
    std::cerr << prefix << msg << "\n" << std::flush;
}
}   // namespace detail

inline void info(const std::string &msg) { detail::write("", msg); }
inline void warn(const std::string &msg) { detail::write("WARN: ", msg); }
inline void error(const std::string &msg) { detail::write("ERROR: ", msg); }

inline void close() {
    log_file.close();
}

}   // namespace logger

#endif
