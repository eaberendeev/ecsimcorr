#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace check {
inline void require(bool value, const std::string& message) {
    if (!value)
        throw std::runtime_error(message);
}
inline void near(double actual, double expected, double tolerance, const std::string& message) {
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
    require(actual == expected, message);
}
// Name the physical/algebraic question in both successful output and failures.
// This wrapper only reports a case; setup, numerical operations and criteria
// remain in the test function in their execution order.
template <typename Test>
void scenario(const std::string& name, Test&& test) {
    std::cout << "[ RUN  ] " << name << '\n';
    try {
        test();
    } catch (const std::exception& error) {
        throw std::runtime_error(name + ": " + error.what());
    }
    std::cout << "[ PASS ] " << name << '\n';
}

template <typename Run>
int run(Run&& test) {
    try {
        test();
        std::cout << "invariants passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
}   // namespace check
