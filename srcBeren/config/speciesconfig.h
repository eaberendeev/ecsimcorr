#pragma once

#include <nlohmann/json.hpp>

#include "particles/species_parameters.h"

SpeciesParameters species_parameters_from_config(const nlohmann::json& config);
