#pragma once

#include <nlohmann/json.hpp>

class Domain;
class BoundaryConditionHandler;

// Build field constraints and particle wall laws after the domain geometry is configured.
// Validate secondary-emission product species after particle initialization.
BoundaryConditionHandler boundary_conditions_from_config(const nlohmann::json& config, const Domain& domain);
