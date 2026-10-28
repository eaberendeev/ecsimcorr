#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace check {
namespace detail {
// Per-process counters: run() resets them before the test, and scenario()
// uses them to reject cases that silently execute no assertions. They are not
// atomic, so check::* must not be called from an OpenMP region.
inline int& check_count() {
    static int count = 0;
    return count;
}
}   // namespace detail

inline void require(bool value, const std::string& message) {
    ++detail::check_count();
    if (!value)
        throw std::runtime_error(message);
}
inline void near(double actual, double expected, double tolerance, const std::string& message) {
    ++detail::check_count();
    if (std::isfinite(actual) && std::isfinite(expected) && std::isfinite(tolerance) && tolerance >= 0 &&
        std::abs(actual - expected) <= tolerance)
        return;
    std::ostringstream details;
    details << std::scientific << std::setprecision(12) << message << ": actual=" << actual << ", expected=" << expected
            << ", tolerance=" << tolerance;
    throw std::runtime_error(details.str());
}
inline void equal(double actual, double expected, double tolerance, const std::string& message) {
    near(actual, expected, tolerance, message);
}
template <typename Actual, typename Expected>
void equal(const Actual& actual, const Expected& expected, const std::string& message) {
    static_assert(!std::is_floating_point_v<std::decay_t<Actual>> &&
                      !std::is_floating_point_v<std::decay_t<Expected>> && !std::is_pointer_v<std::decay_t<Actual>> &&
                      !std::is_pointer_v<std::decay_t<Expected>>,
                  "check::equal with floating point or pointers is ambiguous; use check::near for numbers and an "
                  "explicit string/pointer comparison");
    require(actual == expected, message);
}
// Name the physical/algebraic question in both successful output and failures.
// This wrapper only reports a case; setup, numerical operations and criteria
// remain in the test function in their execution order.
template <typename Test>
void scenario(const std::string& name, Test&& test) {
    std::cout << "[ RUN  ] " << name << '\n';
    const int checks_before = detail::check_count();
    try {
        test();
    } catch (const std::exception& error) {
        throw std::runtime_error(name + ": " + error.what());
    }
    if (detail::check_count() == checks_before)
        throw std::runtime_error(name + ": scenario performed no checks");
    std::cout << "[ PASS ] " << name << '\n';
}

template <typename Run>
int run(Run&& test) {
    detail::check_count() = 0;
    try {
        test();
        if (detail::check_count() == 0)
            throw std::runtime_error("no checks executed");
        std::cout << "invariants passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
}   // namespace check
