#pragma once

#include <limits>
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

// Окно шагов для записи микротаймер-профилей profile_*.json.
// По умолчанию [0, INT_MAX] — профилирование как раньше (вне окна ничего не пишется).
static inline int timerStartStep() {
    const static int res = getenvParsed<int>("TIMER_START_STEP", 0);
    return res;
}

static inline int timerEndStep() {
    const static int res = getenvParsed<int>("TIMER_END_STEP", std::numeric_limits<int>::max());
    return res;
}

// Вариант "path": обрезать траекторию частицы до открытой Z-стенки при депозите
// тока. Z_PATH=0 (по умолчанию) — поведение master; Z_PATH=1 — path.
static inline bool zPathEnabled() {
    const static bool res = getenvParsed<bool>("Z_PATH", false);
    return res;
}

}   // namespace envOptions
