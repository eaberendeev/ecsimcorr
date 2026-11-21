#pragma once
#include <memory>
#include <vector>

#include "grid/domain.h"
#include "utils/containers.h"

struct ElectricFieldConfig {
    virtual ~ElectricFieldConfig() = default;
    // Важно: apply ДОЛЖЕН ДОБАВЛЯТЬ вклад к полю (+=), не затирая существующее
    virtual void apply(Field3d& fieldE, const Domain& domain) const = 0;
};

// Однородное E = (ex,ey,ez)
struct ElectricUniformFieldConfig : ElectricFieldConfig {
    double ex, ey, ez;

    ElectricUniformFieldConfig(double x, double y, double z) : ex(x), ey(y), ez(z) {
    }

    void apply(Field3d& fieldE, const Domain& domain) const override;
};

// Равномерно заряженный цилиндр (вклад в E)
// Plane wave E_c += amplitude_c cos(k . r + phase) at the Yee position of E_c
// (fields/cosine_field.h).
struct ElectricCosineFieldConfig : ElectricFieldConfig {
    Vector3R amplitude, mode;
    double phase = 0.0;

    ElectricCosineFieldConfig(const Vector3R& amplitude, const Vector3R& mode, double phase)
        : amplitude(amplitude), mode(mode), phase(phase) {
    }

    void apply(Field3d& fieldE, const Domain& domain) const override;
};

struct ElectricUniformlyChargedCylinderConfig : ElectricFieldConfig {
    double radius;
    double value;

    ElectricUniformlyChargedCylinderConfig(double r, double v) : radius(r), value(v) {
    }

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
