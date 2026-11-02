#pragma once

#include <nlohmann/json.hpp>

#include "particles/species_state.h"
#include "sources/particle_sources.h"

SpeciesParameters species_parameters_from_config(const nlohmann::json& config);
ParticleSources particle_sources_from_config(const nlohmann::json& config, const SpeciesParameters& parameters,
                                             const Domain& domain);
