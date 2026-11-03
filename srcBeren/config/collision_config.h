#pragma once

#include <nlohmann/json.hpp>

#include "collisions/collision_manager.h"

// Channels retain JSON order. Missing or empty Collisions disables collisions.
CollisionManager collisions_from_config(const nlohmann::json& config, double n0);
