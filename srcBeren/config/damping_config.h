#pragma once

#include <nlohmann/json.hpp>

#include "fields/Damping.h"

// Read once before stepping; restart reconstructs the same layers from the input.
DampingParameters damping_parameters_from_config(const nlohmann::json& config);
