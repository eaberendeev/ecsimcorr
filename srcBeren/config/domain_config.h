#pragma once

#include <nlohmann/json.hpp>

#include "grid/domain.h"

// Configure before constructing particles or field operators that refer to the domain.
// JSON belongs to configuration, so geometry and numerical kernels stay independent of the parser.
Domain domain_from_config(const nlohmann::json& config);

class BoundaryConditionHandler;
BoundaryConditionHandler boundary_conditions_from_config(const nlohmann::json& config, const Domain& domain);
