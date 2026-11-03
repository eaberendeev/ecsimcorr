#include "external_electric_field_config.h"

#include <stdexcept>
#include <utility>

using json = nlohmann::json;

// Парсинг единого блока-конфига (объекта)
static std::unique_ptr<ElectricFieldConfig> parse_electric_object(const json& obj) {
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

std::unique_ptr<ElectricFieldConfig> create_electric_field_config(const json& system_config, const char* key) {
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
            json wrapper;
            wrapper[key] = value;
            if (auto cfg = parse_electric_object(wrapper)) {
                comp->add(std::move(cfg));
            }
        }
        return comp->parts.empty() ? nullptr : std::move(comp);
    }

    return nullptr;
}
