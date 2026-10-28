#include "OutputFactory.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "Read.h"
#include "logger.h"
#include "outputs/EnergySpectrumOutput.h"
#include "outputs/ProbeOutput.h"
#include "outputs/SliceOutputs.h"
#include "outputs/StateOutputs.h"
#include "outputs/TableOutputs.h"

namespace {
[[noreturn]] void invalid_setting(const std::string& name, const std::string& reason) {
    const std::string message = "OutputFactory: '" + name + "' " + reason;
    logger::error(message);
    throw std::runtime_error(message);
}

int read_integer(const nlohmann::json& value, const std::string& name) {
    if (!value.is_number_integer() || value.get<double>() < std::numeric_limits<int>::min() ||
        value.get<double>() > std::numeric_limits<int>::max())
        invalid_setting(name, "must be an integer in the int range");
    return value.get<int>();
}

double read_coordinate(const nlohmann::json& value, const std::string& name) {
    if (!value.is_number() || !std::isfinite(value.get<double>()))
        invalid_setting(name, "must contain finite numeric coordinates");
    return value.get<double>();
}

DiagnosticPlanes parse_planes(const nlohmann::json& config, const Domain& domain) {
    DiagnosticPlanes planes;
    if (!config.contains("planes"))
        return planes;
    const auto& entries = config["planes"];
    if (!entries.is_object())
        invalid_setting("planes", "must be an object with x/y/z arrays");
    for (const auto& entry : entries.items())
        if (entry.key() != "x" && entry.key() != "y" && entry.key() != "z")
            invalid_setting("planes", "has unknown axis '" + entry.key() + "'");
    const std::array<std::vector<double>*, 3> coordinates{&planes.x, &planes.y, &planes.z};
    for (int axis = 0; axis < 3; ++axis) {
        const std::string name(1, "xyz"[axis]);
        if (!entries.contains(name))
            continue;
        if (!entries[name].is_array())
            invalid_setting("planes." + name, "must be an array");
        for (const auto& entry : entries[name]) {
            const double coordinate = read_coordinate(entry, "planes." + name);
            // Match the existing plane index coord/cell_size and keep its
            // out-of-range skip behavior. Check representability before int_value().
            const double index = coordinate / domain.cell_size(axis);
            if (!std::isfinite(index) || index <= std::numeric_limits<int>::min() + 1.0 ||
                index >= std::numeric_limits<int>::max() - 1.0)
                invalid_setting("planes." + name, "cannot be represented as a plane index");
            coordinates[axis]->push_back(coordinate);
        }
    }
    return planes;
}

int get_interval(const nlohmann::json& config, int default_val = 1) {
    const int interval = config.contains("interval") ? read_integer(config["interval"], "interval") : default_val;
    if (interval <= 0)
        invalid_setting("interval", "must be > 0");
    return interval;
}

// For 2D slices, an absent interval means no instantaneous output.
int get_slice_interval(const nlohmann::json& config) {
    return config.contains("interval") ? get_interval(config) : 0;
}

int get_average_steps(const nlohmann::json& config) {
    const int steps = config.contains("average_steps") ? read_integer(config["average_steps"], "average_steps") : 0;
    if (steps < 0)
        invalid_setting("average_steps", "must be >= 0");
    return steps;
}

std::vector<std::string> resolve_species_names(const Species& all_species, const nlohmann::json& config) {
    std::vector<std::string> names;
    if (!config.contains("species") || (config["species"].is_string() && config["species"] == "all")) {
        for (const auto& kv : all_species) names.push_back(kv.first);
        return names;
    }
    const auto& entries = config["species"];
    if (!entries.is_array())
        invalid_setting("species", "must be 'all' or an array of species names");
    for (const auto& entry : entries) {
        if (!entry.is_string())
            invalid_setting("species", "must contain only species names");
        const auto name = entry.get<std::string>();
        if (!all_species.contains(name))
            invalid_setting("species", "contains unknown species '" + name + "'");
        names.push_back(name);
    }
    return names;
}

SliceOutputSettings parse_slice_settings(const nlohmann::json& config, const Domain& domain) {
    const int interval = get_slice_interval(config);
    const int average_steps = get_average_steps(config);
    return {parse_planes(config, domain), interval, average_steps};
}

std::vector<std::string> parse_field_names(const nlohmann::json& config) {
    if (!config.contains("fields"))
        return {"E", "B"};
    const auto& entries = config["fields"];
    if (!entries.is_array())
        invalid_setting("fields", "must be an array of field names");
    const bool slices = config["type"] == "fields_2d";
    std::vector<std::string> names;
    for (const auto& entry : entries) {
        if (!entry.is_string())
            invalid_setting("fields", "must contain only field names");
        const auto name = entry.get<std::string>();
        if (name != "E" && name != "B" && name != "En" && name != "Bn" && name != "E_ex" &&
            !(slices && name == "B_init"))
            invalid_setting("fields",
                            "contains unsupported field '" + name + "' for " + config["type"].get<std::string>());
        names.push_back(name);
    }
    return names;
}

std::vector<int> parse_timesteps(const nlohmann::json& config) {
    std::vector<int> steps;
    if (!config.contains("timesteps"))
        return steps;
    if (!config["timesteps"].is_array())
        invalid_setting("timesteps", "must be an array");
    for (const auto& entry : config["timesteps"]) {
        const int step = read_integer(entry, "timesteps");
        if (step < 0)
            invalid_setting("timesteps", "must contain nonnegative steps");
        steps.push_back(step);
    }
    return steps;
}

std::vector<Vector3R> parse_probe_points(const nlohmann::json& config, const Domain& domain) {
    std::vector<Vector3R> points;
    if (!config.contains("points"))
        return points;
    if (!config["points"].is_array())
        invalid_setting("points", "must be an array of [x,y,z] points");
    for (const auto& entry : config["points"]) {
        if (!entry.is_array() || entry.size() != 3)
            invalid_setting("points", "must contain [x,y,z] points");
        const Vector3R point(read_coordinate(entry[0], "points"), read_coordinate(entry[1], "points"),
                             read_coordinate(entry[2], "points"));
        const auto normalized = domain.to_cell_coordinates(point);
        for (int axis = 0; axis < 3; ++axis) {
            // Linear scalar/E/B gather uses neighboring nodes at integer and
            // half-integer positions. Validate their stencil, including ghosts;
            // probes at a physical wall remain valid.
            const double node = normalized[axis] + GHOST_CELLS;
            if (!std::isfinite(node) || node < 0.5 || node >= domain.size()[axis] - 1)
                invalid_setting("points", "has a probe outside the interpolation stencil");
        }
        points.push_back(point);
    }
    return points;
}

}   // namespace

// ============================================================
//  OutputFactory
// ============================================================
std::vector<std::unique_ptr<IDiagnosticOutput>> OutputFactory::create(const nlohmann::json& diag_config,
                                                                      const Domain& domain, const Species& species,
                                                                      const nlohmann::json& system_config,
                                                                      const Field3d& fieldEn, const Field3d& fieldBn,
                                                                      const Field3d& fieldE_external,
                                                                      const Field3d& fieldBInit) {
    if (!diag_config.is_object())
        invalid_setting("diagnostics", "must be an object");
    if (!diag_config.contains("outputs"))
        return {};
    if (!diag_config["outputs"].is_array()) {
        throw std::runtime_error(
            "OutputFactory: missing 'outputs' array in diagnostics config. "
            "Update your gen_config.py to use the new format.");
    }
    return create_from_list(diag_config["outputs"], domain, species, system_config, fieldEn, fieldBn, fieldE_external,
                            fieldBInit);
}

std::vector<std::unique_ptr<IDiagnosticOutput>> OutputFactory::create_from_list(
    const nlohmann::json& outputs, const Domain& domain, const Species& species, const nlohmann::json& system_config,
    const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external, const Field3d& fieldBInit) {
    std::vector<std::unique_ptr<IDiagnosticOutput>> result;

    for (const auto& cfg : outputs) {
        if (!cfg.is_object() || !cfg.contains("type") || !cfg["type"].is_string())
            invalid_setting("type", "must be specified as a string in every output object");
        const std::string type = cfg["type"].get<std::string>();

        if (type == "energy_balance") {
            result.push_back(
                std::make_unique<EnergyBalanceOutput>(get_checked<double>(system_config, "Dt"), get_interval(cfg)));
        } else if (type == "boundary_stats") {
            result.push_back(
                std::make_unique<BoundaryOutput>(get_checked<double>(system_config, "Dt"), get_interval(cfg)));
        } else if (type == "console_summary") {
            result.push_back(std::make_unique<ConsoleOutput>(species, get_checked<double>(system_config, "Dt")));
        } else if (type == "fields_2d") {
            result.push_back(std::make_unique<Fields2DOutput>(domain, fieldEn, fieldBn, fieldE_external, fieldBInit,
                                                              parse_slice_settings(cfg, domain),
                                                              parse_field_names(cfg)));
        } else if (type == "density_2d") {
            result.push_back(std::make_unique<Density2DOutput>(species, domain, parse_slice_settings(cfg, domain),
                                                               resolve_species_names(species, cfg)));
        } else if (type == "current_2d") {
            result.push_back(std::make_unique<Current2DOutput>(species, domain, parse_slice_settings(cfg, domain),
                                                               resolve_species_names(species, cfg)));
        } else if (type == "pressure_2d") {
            result.push_back(std::make_unique<Pressure2DOutput>(species, domain, parse_slice_settings(cfg, domain),
                                                                resolve_species_names(species, cfg)));
        } else if (type == "energy_spectrum") {
            result.push_back(std::make_unique<EnergySpectrumOutput>(species, resolve_species_names(species, cfg),
                                                                    get_interval(cfg)));
        } else if (type == "recovery") {
            // Recovery is optional: absent interval disables it, but an
            // explicitly configured interval must be positive.
            if (cfg.contains("interval")) {
                const int interval = get_interval(cfg);
                result.push_back(std::make_unique<RecoveryOutput>(fieldEn, fieldBn, species, interval));
            }
        } else if (type == "fields_3d") {
            result.push_back(std::make_unique<Fields3DOutput>(fieldEn, fieldBn, fieldE_external, fieldBInit,
                                                              parse_field_names(cfg), parse_timesteps(cfg)));
        } else if (type == "charge_density_2d") {
            result.push_back(
                std::make_unique<ChargeDensity2DOutput>(species, domain, parse_planes(cfg, domain), get_interval(cfg)));
        } else if (type == "total_current_2d") {
            result.push_back(
                std::make_unique<TotalCurrent2DOutput>(species, domain, parse_planes(cfg, domain), get_interval(cfg)));
        } else if (type == "full3d_all") {
            result.push_back(std::make_unique<Full3DAllOutput>(fieldEn, fieldBn, fieldE_external, fieldBInit, species,
                                                               parse_timesteps(cfg)));
        } else if (type == "probes") {
            result.push_back(std::make_unique<ProbeOutput>(domain, fieldEn, fieldBn, fieldE_external, fieldBInit,
                                                           species, resolve_species_names(species, cfg),
                                                           parse_probe_points(cfg, domain), get_interval(cfg)));
        } else {
            const std::string msg = "OutputFactory: unknown output type '" + type + "'";
            logger::error(msg);
            throw std::runtime_error(msg);
        }
    }
    return result;
}
