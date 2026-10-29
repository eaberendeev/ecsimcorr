#pragma once

#include <nlohmann/json.hpp>

#include "grid/Domain.h"

// Configure before constructing particles or field operators that refer to the domain.
void configure_domain(Domain& domain, const nlohmann::json& config);
