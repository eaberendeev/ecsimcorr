#include "speciesconfig.h"

#include <cmath>
#include <stdexcept>

SpeciesParameters species_parameters_from_config(const nlohmann::json& config) {
    SpeciesParameters parameters{config.value("Name", ""),
                                 config.value("Charge", 1.0),
                                 config.value("Mass", 1.0),
                                 config.value("Density", 1.0),
                                 config.at("NumPartPerCell").get<int>(),
                                 0};
    if (parameters.particles_per_cell <= 0 || !std::isfinite(parameters.mass) || parameters.mass <= 0 ||
        !std::isfinite(parameters.density) || parameters.density < 0 || !std::isfinite(parameters.charge))
        throw std::invalid_argument(
            "species requires positive mass/NumPartPerCell and finite charge/nonnegative density");
    parameters.macro_weight = parameters.density / parameters.particles_per_cell;
    return parameters;
}
