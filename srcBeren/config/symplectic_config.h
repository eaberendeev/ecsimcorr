#pragma once

#include <nlohmann/json.hpp>

#include "core/vector3.h"

Vector3R initial_symplectic_vector(const nlohmann::json& config, const char* name);
void validate_symplectic_config(const nlohmann::json& system, const nlohmann::json& particles);
