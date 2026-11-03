#pragma once

#include <nlohmann/json.hpp>

#include "collisions/collision_manager.h"

// Modern channels retain JSON order and suppress the legacy Collider fallback.
CollisionManager collisions_from_config(const nlohmann::json& config, double n0);
