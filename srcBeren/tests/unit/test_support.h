#pragma once

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "ParticlesArray.h"

namespace check {
inline void require(bool value, const std::string& message) {
    if (!value)
        throw std::runtime_error(message);
}
inline void near(double actual, double expected, double tolerance, const std::string& message) {
    require(std::isfinite(actual) && std::isfinite(expected) && std::abs(actual - expected) <= tolerance,
            message + ": actual=" + std::to_string(actual) + ", expected=" + std::to_string(expected));
}
inline void vector(const Vector3R& actual, const Vector3R& expected, double tolerance, const std::string& message) {
    for (int d = 0; d < 3; ++d) near(actual[d], expected[d], tolerance, message);
}
inline nlohmann::json species(const std::string& name, double charge = -1, double mass = 1, double density = 1,
                              int npc = 8) {
    return {{"Name", name},       {"Charge", charge},      {"Mass", mass},
            {"Density", density}, {"NumPartPerCell", npc}, {"distribution", nlohmann::json::array()}};
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
