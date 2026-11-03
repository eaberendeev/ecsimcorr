#include "external_magnetic_field_config.h"

#include <stdexcept>
#include <utility>

using json = nlohmann::json;

// Парсинг единичного объекта-конфига
static std::unique_ptr<MagneticFieldConfig> parse_magnetic_object(const json& obj) {
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

std::unique_ptr<MagneticFieldConfig> create_magnetic_field_config(const json& system_config, const char* key) {
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
            json wrapper;
            wrapper[key] = value;
            if (auto cfg = parse_magnetic_object(wrapper)) {
                comp->add(std::move(cfg));
            }
        }
        return comp->parts.empty() ? nullptr : std::move(comp);
    }

    return nullptr;
}
