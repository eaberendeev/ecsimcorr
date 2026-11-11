#include "field_config.h"

#include <stdexcept>
#include <string>
#include <utility>

namespace {
struct CosineParameters {
    Vector3R amplitude, mode;
    double phase = 0.0;
};

// {"amplitude": [x, y, z], "mode": [x, y, z], "phase": p} of a cosine_field.
CosineParameters parse_cosine(const nlohmann::json& u) {
    const auto vector3 = [&](const char* name) {
        if (!u.contains(name) || !u.at(name).is_array() || u.at(name).size() != 3)
            throw std::invalid_argument(std::string("cosine_field requires ") + name + " = [x, y, z]");
        return Vector3R(u.at(name)[0].get<double>(), u.at(name)[1].get<double>(), u.at(name)[2].get<double>());
    };
    return {vector3("amplitude"), vector3("mode"), u.value("phase", 0.0)};
}
}   // namespace

// Парсинг единого блока-конфига (объекта)
static std::unique_ptr<ElectricFieldConfig> parse_electric_object(const nlohmann::json& obj) {
    if (obj.contains("cosine_field")) {
        const CosineParameters p = parse_cosine(obj.at("cosine_field"));
        return std::make_unique<ElectricCosineFieldConfig>(p.amplitude, p.mode, p.phase);
    }

    if (obj.contains("uniform_field")) {
        const auto& u = obj.at("uniform_field");
        double ex = 0, ey = 0, ez = 0;
        if (u.contains("value") && u.at("value").is_array() && u.at("value").size() == 3) {
            ex = u.at("value")[0].get<double>();
            ey = u.at("value")[1].get<double>();
            ez = u.at("value")[2].get<double>();
        }
        return std::make_unique<ElectricUniformFieldConfig>(ex, ey, ez);
    }

    if (obj.contains("uniformly_charged_cylinder")) {
        const auto& c = obj.at("uniformly_charged_cylinder");
        double radius = c.at("radius").get<double>();
        double value = c.at("value").get<double>();
        return std::make_unique<ElectricUniformlyChargedCylinderConfig>(radius, value);
    }

    return nullptr;
}

std::unique_ptr<ElectricFieldConfig> create_electric_field_config(const nlohmann::json& system_config,
                                                                  const char* key) {
    if (!system_config.contains(key))
        return nullptr;
    const auto& ext = system_config.at(key);

    // Вариант: массив — явный порядок суммирования
    if (ext.is_array()) {
        auto comp = std::make_unique<ElectricCompositeFieldConfig>();
        for (const auto& item : ext) {
            if (!item.is_object())
                continue;
            if (auto cfg = parse_electric_object(item))
                comp->add(std::move(cfg));
        }
        return comp->parts.empty() ? nullptr : std::move(comp);
    }

    if (ext.is_object()) {
        auto comp = std::make_unique<ElectricCompositeFieldConfig>();
        // Перебираем все ключи объекта, каждый обрабатываем как отдельную
        // конфигурацию
        for (auto& [key, value] : ext.items()) {
            nlohmann::json wrapper;
            wrapper[key] = value;
            if (auto cfg = parse_electric_object(wrapper)) {
                comp->add(std::move(cfg));
            }
        }
        return comp->parts.empty() ? nullptr : std::move(comp);
    }

    return nullptr;
}

// Парсинг единичного объекта-конфига
static std::unique_ptr<MagneticFieldConfig> parse_magnetic_object(const nlohmann::json& obj) {
    if (obj.contains("uniform_field")) {
        const auto& u = obj.at("uniform_field");
        double bx = 0, by = 0, bz = 0;
        if (u.contains("value") && u.at("value").is_array() && u.at("value").size() == 3) {
            bx = u.at("value")[0].get<double>();
            by = u.at("value")[1].get<double>();
            bz = u.at("value")[2].get<double>();
        }
        return std::make_unique<MagneticUniformFieldConfig>(bx, by, bz);
    }

    if (obj.contains("cosine_field")) {
        const CosineParameters p = parse_cosine(obj.at("cosine_field"));
        return std::make_unique<MagneticCosineFieldConfig>(p.amplitude, p.mode, p.phase);
    }

    if (obj.contains("coils")) {
        const auto& arr = obj.at("coils");
        if (arr.is_array()) {
            std::vector<Coil> coils;
            for (const nlohmann::json& coil : arr) {
                if (!coil.contains("z") || !coil.contains("R") || !coil.contains("I")) {
                    throw std::invalid_argument("Invalid coil configuration. Missing z, R, or I parameters.");
                }
                coils.emplace_back(coil.at("z").get<double>(), coil.at("R").get<double>(), coil.at("I").get<double>());
            }
            return std::make_unique<MagneticCoilsFieldConfig>(std::move(coils));
        }
    }

    return nullptr;
}

std::unique_ptr<MagneticFieldConfig> create_magnetic_field_config(const nlohmann::json& system_config,
                                                                  const char* key) {
    if (!system_config.contains(key))
        return nullptr;
    const auto& ext = system_config.at(key);

    if (ext.is_array()) {
        auto comp = std::make_unique<MagneticCompositeFieldConfig>();
        for (const auto& item : ext) {
            if (!item.is_object())
                continue;
            if (auto cfg = parse_magnetic_object(item))
                comp->add(std::move(cfg));
        }
        return comp->parts.empty() ? nullptr : std::move(comp);
    }

    if (ext.is_object()) {
        auto comp = std::make_unique<MagneticCompositeFieldConfig>();
        // Перебираем все ключи объекта, каждый обрабатываем как отдельную
        // конфигурацию
        for (auto& [key, value] : ext.items()) {
            nlohmann::json wrapper;
            wrapper[key] = value;
            if (auto cfg = parse_magnetic_object(wrapper)) {
                comp->add(std::move(cfg));
            }
        }
        return comp->parts.empty() ? nullptr : std::move(comp);
    }

    return nullptr;
}

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
