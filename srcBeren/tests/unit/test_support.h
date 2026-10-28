#pragma once

#include "ParticlesArray.h"
#include "test_case.h"

namespace check {
inline void vector(const Vector3R& actual, const Vector3R& expected, double tolerance, const std::string& message) {
    for (int d = 0; d < 3; ++d) near(actual[d], expected[d], tolerance, message);
}
inline nlohmann::json species(const std::string& name, double charge = -1, double mass = 1, double density = 1,
                              int npc = 8) {
    return {{"Name", name},       {"Charge", charge},      {"Mass", mass},
            {"Density", density}, {"NumPartPerCell", npc}, {"distribution", nlohmann::json::array()}};
}
}   // namespace check
