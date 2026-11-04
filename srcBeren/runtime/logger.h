#pragma once

#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// Logging is a serial-only facility: init/close and all messages are emitted
// from the orchestration level, never from inside an OpenMP parallel region.
// Locking here would serialise the hot loop for messages that must not happen
// there anyway, so the writers stay unsynchronised on purpose.
namespace logger {

inline std::ofstream log_file;
inline int current_timestep = 0;

inline void init(const std::string& path = "beren3d.log") {
    log_file.open(path);
    log_file << "=== beren3d log ===\n\n";
}

inline void set_timestep(int ts) {
    current_timestep = ts;
}

namespace detail {
inline void write(const std::string& prefix, const std::string& msg, bool flush) {
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

inline void info(const std::string& msg) {
    detail::write("", msg, /*flush=*/false);
}
inline void warn(const std::string& msg) {
    detail::write("WARN: ", msg, /*flush=*/true);
}
inline void error(const std::string& msg) {
    detail::write("ERROR: ", msg, /*flush=*/true);
}

inline void close() {
    log_file.flush();
    log_file.close();
}

}   // namespace logger

#endif
