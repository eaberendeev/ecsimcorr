#pragma once

#include <string>

#include "util.h"

namespace envOptions {
template <typename T>
static inline T getenvParsed(const char* name, T def);

template <>
inline int getenvParsed(const char* name, int def) {
    const char* value = getenv(name);
    if (value == nullptr) {
        return def;
    }

    return std::stoi(value);
}

template <>
inline int64_t getenvParsed(const char* name, int64_t def) {
    const char* value = getenv(name);
    if (value == nullptr) {
        return def;
    }

    return std::stoll(value);
}

template <>
inline bool getenvParsed(const char* name, bool def) {
    const char* value = getenv(name);
    if (value == nullptr) {
        return def;
    }

    return static_cast<bool>(std::stoi(value));
}

static inline int validationPeriodicity() {
    const static int res = getenvParsed<int>("VALIDATION_PERIODICITY", 10);
    return res;
}

static inline int timeredSimSteps() {
    const static int res = getenvParsed<int>("TIMERED_SIM_STEPS", 20);
    return res;
}

static inline bool useMixedPrecision() {
    const static bool res = getenvParsed<bool>("USE_MIXED_PRECISION", true);
    return res;
}
static inline int64_t timerThresholdNs() {
    const static int64_t res = getenvParsed<int64_t>("TIMER_THRESHOLD_NS", -1);
    return res;
}

}   // namespace envOptions
