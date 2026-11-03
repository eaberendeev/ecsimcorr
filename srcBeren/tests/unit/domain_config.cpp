#include "config/domain_config.h"

#include <stdexcept>
#include <string>

#include "test_case.h"

void domain_configuration_preserves_geometry() {
    nlohmann::json config = {{"Dx", 1}, {"Dy", 1}, {"Dz", 1}, {"NumCellsX", 4}, {"NumCellsY", 4}, {"NumCellsZ", 4}};
    Domain domain = domain_from_config(config);
    check::require(!domain.geom.use_cylinder, "box config enabled a cylinder");
    for (int axis = 0; axis < 3; ++axis) check::near(domain.geom.box_max[axis], 4, 0, "configured box dimensions");
    config["CylinderDomain"] = {{"radius", 0.5}, {"center", {1, 2, 7}}};
    domain = domain_from_config(config);
    check::require(domain.geom.use_cylinder, "cylinder config was ignored");
    // Existing input uses only transverse center coordinates; z stays zero.
    const Vector3R expected_center(1, 2, 0);
    for (int axis = 0; axis < 3; ++axis)
        check::near(domain.geom.cyl_center[axis], expected_center[axis], 0, "configured cylinder center");
    check::near(domain.geom.cyl_radius, 0.5, 0, "configured cylinder radius");
    config["CylinderDomain"].erase("center");
    domain = domain_from_config(config);
    check::require(!domain.geom.use_cylinder, "legacy radius-only config must leave the cylinder disabled");
    config.erase("Dx");
    bool rejected = false;
    try {
        domain = domain_from_config(config);
    } catch (const std::runtime_error& error) {
        rejected = std::string(error.what()).find("Dx") != std::string::npos;
    }
    check::require(rejected, "domain configuration must still identify a missing Dx");
}

int main() {
    return check::run([] {
        check::scenario("Domain configuration preserves box/cylinder input and missing-key errors",
                        domain_configuration_preserves_geometry);
    });
}
