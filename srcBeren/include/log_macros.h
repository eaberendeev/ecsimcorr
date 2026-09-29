#pragma once

#include <omp.h>

#include <sstream>

#include "logger.h"
#include "timer.h"

inline bool g_verbose_step = false;

#define LOG_STEP(x)                                        \
    do {                                                   \
        if (g_verbose_step && omp_get_thread_num() == 0) { \
            assert(!omp_in_parallel());                    \
            timer::commonTimer __timerVerbose("LOG_STEP"); \
            std::ostringstream _ls;                        \
            _ls << x;                                      \
            logger::info(_ls.str());                       \
        }                                                  \
    } while (0)
