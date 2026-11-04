#pragma once

#include <sstream>

#include "runtime/logger.h"
#include "runtime/profiler.h"

inline bool g_verbose_step = false;

#define LOG_STEP(x)                                        \
    do {                                                   \
        if (g_verbose_step) {                              \
            timer::commonTimer __timerVerbose("LOG_STEP"); \
            std::ostringstream _ls;                        \
            _ls << x;                                      \
            logger::info(_ls.str());                       \
        }                                                  \
    } while (0)
