#pragma once

#include <sstream>

#include "logger.h"
#include "timer.h"

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
