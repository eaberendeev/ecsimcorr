#include "symplectic_config.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

Vector3R initial_symplectic_vector(const nlohmann::json& config, const char* name) {
    const std::vector<double> values = config.value(name, std::vector<double>{0, 0, 0});
    if (values.size() != 3)
        throw std::invalid_argument(std::string(name) + " requires three components");
    for (double x : values)
        if (!std::isfinite(x))
            throw std::invalid_argument(std::string(name) + " must be finite");
    return {values[0], values[1], values[2]};
}

void validate_symplectic_config(const nlohmann::json& system, const nlohmann::json& particles) {
    // Validate the Hamiltonian model before allocating the mesh or loading particles.
    const double dt = system.at("Dt").get<double>();
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("symplectic requires finite positive Dt");

    double inverse_h2 = 0;
    for (const char* key : {"Dx", "Dy", "Dz"}) {
        const double h = system.at(key).get<double>();
        if (!std::isfinite(h) || h <= 0)
            throw std::invalid_argument("symplectic requires positive cell sizes");
        inverse_h2 += 1 / (h * h);
    }
    if (dt * std::sqrt(inverse_h2) >= 1)
        throw std::invalid_argument("symplectic explicit Maxwell CFL requires Dt*sqrt(sum(1/h^2)) < 1");

    for (const char* key : {"NumCellsX", "NumCellsY", "NumCellsZ"})
        if (system.at(key).get<int>() < 1)
            throw std::invalid_argument("symplectic requires nonempty grid");

    if (system.contains("CylinderDomain") || system.value("DampingType", "None") != "None" ||
        system.contains("Collider") || (system.contains("Collisions") && !system.at("Collisions").empty()))
        throw std::invalid_argument("symplectic supports collisionless periodic boxes without damping");

    for (const char* key : {"ExternalFieldE", "ExternalFieldB"})
        if (system.contains(key) && !system.at(key).empty())
            throw std::invalid_argument(
                "symplectic: use Symplectic.initial_electric/initial_magnetic, not external fields");

    for (const nlohmann::json& boundary : system.at("Boundary_conditions"))
        if (boundary.size() != 1 || !boundary.contains("periodic"))
            throw std::invalid_argument("symplectic supports only periodic boundary conditions");

    const nlohmann::json config = system.value("Symplectic", nlohmann::json::object());
    for (nlohmann::json::const_iterator it = config.begin(); it != config.end(); ++it)
        if (it.key() != "initial_electric" && it.key() != "initial_magnetic")
            throw std::invalid_argument("unknown Symplectic option: " + it.key());
    initial_symplectic_vector(config, "initial_electric");
    initial_symplectic_vector(config, "initial_magnetic");

    for (const nlohmann::json& sp : particles.at("particles"))
        for (const nlohmann::json& dist : sp.value("distribution", nlohmann::json::array()))
            if (dist.value("type", "") != "initial")
                throw std::invalid_argument("symplectic does not support injection");
}
