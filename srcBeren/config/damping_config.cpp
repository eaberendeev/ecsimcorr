#include "damping_config.h"

#include <stdexcept>
#include <string>

DampingParameters damping_parameters_from_config(const nlohmann::json& config) {
    DampingParameters parameters;
    const std::string type = config.value("DampingType", "None");
    if (type == "None")
        return parameters;
    if (type == "CircleXY")
        parameters.type = DampingType::CircleXY;
    else if (type == "Rectangle")
        parameters.type = DampingType::Rectangle;
    else
        throw std::invalid_argument("DampingType is not defined: " + type);

    parameters.coefficient = config.value("DampingCoeff", 0.8);
    if (parameters.type == DampingType::CircleXY) {
        const auto widths = config.value("DampCellsX_glob", nlohmann::json::array({0, 0}));
        if (widths.is_array() && !widths.empty())
            parameters.x_cells[0] = widths[0].get<int>();
        if (widths.is_array() && widths.size() > 1)
            parameters.x_cells[1] = widths[1].get<int>();
    } else {
        parameters.x_cells[0] = config.at("DampCellsX_glob")[0].get<int>();
        parameters.y_cells[0] = config.at("DampCellsY_glob")[0].get<int>();
        parameters.z_cells[0] = config.at("DampCellsZ_glob")[0].get<int>();
        parameters.x_cells[1] = config.at("DampCellsX_glob")[1].get<int>();
        parameters.y_cells[1] = config.at("DampCellsY_glob")[1].get<int>();
        parameters.z_cells[1] = config.at("DampCellsZ_glob")[1].get<int>();
    }
    return parameters;
}
