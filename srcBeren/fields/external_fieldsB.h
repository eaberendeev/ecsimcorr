#pragma once
#include <memory>
#include <vector>
#include "containers.h"
#include "World.h"
#include "nlohmann/json.hpp"

struct MagneticFieldConfig {
    virtual ~MagneticFieldConfig() = default;
    // Принцип: apply ДОЛЖЕН ДОБАВЛЯТЬ вклад (+=)
    virtual void apply(Field3d& fieldB, const Domain& domain) const = 0;
};

// Однородное B = (bx,by,bz)
struct MagneticUniformFieldConfig : MagneticFieldConfig {
    double bx, by, bz;
    MagneticUniformFieldConfig(double x, double y, double z)
        : bx(x), by(y), bz(z) {}
    void apply(Field3d& fieldB, const Domain& domain) const override;
};

// Катушки (суммирует вклад катушек)
struct MagneticCoilsFieldConfig : MagneticFieldConfig {
    nlohmann::json coils_json;   // ожидается массив объектов {z,R,I}
    explicit MagneticCoilsFieldConfig(nlohmann::json coils)
        : coils_json(std::move(coils)) {}
    void apply(Field3d& fieldB, const Domain& domain) const override;
};

// Bphi через Jz
struct BphiConfig final : MagneticFieldConfig {
    double Jz, radius, startz, endz;
    BphiConfig(double j, double r, double s, double e)
        : Jz(j), radius(r), startz(s), endz(e) {}
    void apply(Field3d& fieldB, const Domain& domain) const override;
};

// Композит — сумма нескольких магнитных конфигураций
struct MagneticCompositeFieldConfig : MagneticFieldConfig {
    std::vector<std::unique_ptr<MagneticFieldConfig>> parts;
    void add(std::unique_ptr<MagneticFieldConfig> p) {
        parts.emplace_back(std::move(p));
    }
    void apply(Field3d& fieldB, const Domain& domain) const override;
};

void set_Bphi(Field3d& fieldB, double Jz, double radius, double startz,
              double endz, const Domain& domain);
// Фабрика для ExternalFieldB
// Поддерживается:
//  - { "uniform_field": { "value": [bx,by,bz] } }
//  - { "coils": [ {z,R,I}, ... ] }
//  - { "Bphi": [ {"Jz" : jz}, {"radius" : radius} ] }
//  - массив [{...}, {...}] для явного порядка
std::unique_ptr<MagneticFieldConfig> create_magnetic_field_config(
    const nlohmann::json& system_config, const char* key = "ExternalFieldB");