#pragma once

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>

#include "util.h"

namespace envOptions {

namespace detail {

// Parses an integer string, accepting leading and trailing whitespace but
// rejecting trailing garbage and values outside the range of Int. Throws
// std::invalid_argument with the variable name and offending value on failure.
template <typename Int>
inline Int parseIntegerEnv(const char* name, const char* rawValue) {
    const std::string raw(rawValue);
    const char* begin = raw.data();
    const char* end = begin + raw.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(*begin))) ++begin;
    while (end > begin && std::isspace(static_cast<unsigned char>(end[-1]))) --end;

    // std::stoll/std::stoul accepted a leading '+', while std::from_chars does not.
    if (end - begin >= 2 && begin[0] == '+' && begin[1] >= '0' && begin[1] <= '9')
        ++begin;

    Int parsed = 0;
    const auto result = std::from_chars(begin, end, parsed);
    if (result.ec == std::errc::invalid_argument) {
        throw std::invalid_argument(std::string("Environment variable ") + name + " must be an integer, got \"" +
                                    rawValue + "\"");
    }
    if (result.ec == std::errc::result_out_of_range) {
        throw std::invalid_argument(std::string("Environment variable ") + name +
                                    " is out of range for its integer type, got \"" + rawValue + "\"");
    }
    if (result.ptr != end) {
        throw std::invalid_argument(std::string("Environment variable ") + name +
                                    " must be an integer without trailing characters, got \"" + rawValue + "\"");
    }
    return parsed;
}

}   // namespace detail

template <typename T>
static inline T getenvParsed(const char* name, T def);

template <>
inline int getenvParsed(const char* name, int def) {
    const char* value = getenv(name);
    if (value == nullptr) {
        return def;
    }

    return detail::parseIntegerEnv<int>(name, value);
}

template <>
inline int64_t getenvParsed(const char* name, int64_t def) {
    const char* value = getenv(name);
    if (value == nullptr) {
        return def;
    }

    return detail::parseIntegerEnv<int64_t>(name, value);
}

template <>
inline bool getenvParsed(const char* name, bool def) {
    const char* value = getenv(name);
    if (value == nullptr) {
        return def;
    }

    std::string lowered(value);
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const auto isSpace = [](unsigned char c) { return std::isspace(c) != 0; };
    lowered.erase(lowered.begin(), std::find_if_not(lowered.begin(), lowered.end(), isSpace));
    lowered.erase(std::find_if_not(lowered.rbegin(), lowered.rend(), isSpace).base(), lowered.end());
    if (lowered == "1" || lowered == "true" || lowered == "yes") {
        return true;
    }
    if (lowered == "0" || lowered == "false" || lowered == "no") {
        return false;
    }
    throw std::invalid_argument(std::string("Environment variable ") + name +
                                " must be a boolean (0/1/true/false/yes/no), got \"" + value + "\"");
}

static inline int validationPeriodicity() {
    const static int res = getenvParsed<int>("VALIDATION_PERIODICITY", 10);
    if (res <= 0)
        throw std::invalid_argument("VALIDATION_PERIODICITY must be positive");
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

// Deterministic base seed for Coulomb collision streams (default matches the
// historical hard-coded value).  Overridable for seed-ensemble studies.
static inline int64_t collisionSeed() {
    const static int64_t res = getenvParsed<int64_t>("COLLISION_SEED", 42);
    return res;
}

}   // namespace envOptions
