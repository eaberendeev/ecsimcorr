#include "domainconfig.h"

#include <iostream>
#include <nlohmann/json.hpp>

#include "jsonsettings.h"

Domain domain_from_config(const nlohmann::json& config) {
    Domain domain;
    const Vector3R cell_size = Vector3R(get_checked<double>(config, "Dx"), get_checked<double>(config, "Dy"),
                                        get_checked<double>(config, "Dz"));

    const Vector3I num_cells = Vector3I(get_checked<int>(config, "NumCellsX"), get_checked<int>(config, "NumCellsY"),
                                        get_checked<int>(config, "NumCellsZ"));

    domain.init(num_cells, cell_size);
    if (config.contains("CylinderDomain") && config["CylinderDomain"].is_object()) {
        const nlohmann::json& cyl = config["CylinderDomain"];
        const double radius = cyl["radius"].get<double>();
        if (cyl.contains("center")) {
            const Vector3R center = Vector3R(cyl["center"][0], cyl["center"][1], 0);
            std::cout << "cyl: radius " << radius << ". center: " << center.x() << ", " << center.y() << "\n";
            domain.set_cylinder(center, radius);
        }
    }
    return domain;
}
