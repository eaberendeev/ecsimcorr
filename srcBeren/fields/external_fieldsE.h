#pragma once
#include <memory>
#include <vector>

#include "World.h"
#include "containers.h"
#include "nlohmann/json.hpp"

struct ElectricFieldConfig {
    virtual ~ElectricFieldConfig() = default;
    // Важно: apply ДОЛЖЕН ДОБАВЛЯТЬ вклад к полю (+=), не затирая существующее
    virtual void apply(Field3d& fieldE, const Domain& domain) const = 0;
};

// Однородное E = (ex,ey,ez)
struct ElectricUniformFieldConfig : ElectricFieldConfig {
    double ex, ey, ez;
    ElectricUniformFieldConfig(double x, double y, double z)
        : ex(x), ey(y), ez(z) {}
    void apply(Field3d& fieldE, const Domain& domain) const override;
};

// Равномерно заряженный цилиндр (вклад в E)
struct ElectricUniformlyChargedCylinderConfig : ElectricFieldConfig {
    double radius;
    double value;
    ElectricUniformlyChargedCylinderConfig(double r, double v)
        : radius(r), value(v) {}
    void apply(Field3d& fieldE, const Domain& domain) const override;
};

// Композит — сумма нескольких электрических конфигураций
struct ElectricCompositeFieldConfig : ElectricFieldConfig {
    std::vector<std::unique_ptr<ElectricFieldConfig>> parts;
    void add(std::unique_ptr<ElectricFieldConfig> p) {
        parts.emplace_back(std::move(p));
    }
    void apply(Field3d& fieldE, const Domain& domain) const override;
};

// Фабрика для ExternalFieldE: принимает либо объект, либо массив объектов
// Поддерживается:
//  - { "uniform_field": { "value": [ex,ey,ez] } }
//  - { "uniformly_charged_cylinder": { "radius": R, "value": V } }
//  - или массив [{...}, {...}]
std::unique_ptr<ElectricFieldConfig> create_electric_field_config(
    const nlohmann::json& system_config, const char* key = "ExternalFieldE");