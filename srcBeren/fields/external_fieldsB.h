#pragma once
#include <memory>
#include <vector>

#include "Coil.h"
#include "grid/domain.h"
#include "utils/containers.h"

struct MagneticFieldConfig {
    virtual ~MagneticFieldConfig() = default;
    // Принцип: apply ДОЛЖЕН ДОБАВЛЯТЬ вклад (+=)
    virtual void apply(Field3d& fieldB, const Domain& domain) const = 0;
};

// Однородное B = (bx,by,bz)
struct MagneticUniformFieldConfig : MagneticFieldConfig {
    double bx, by, bz;

    MagneticUniformFieldConfig(double x, double y, double z) : bx(x), by(y), bz(z) {
    }

    void apply(Field3d& fieldB, const Domain& domain) const override;
};

// Катушки (суммирует вклад катушек)
struct MagneticCoilsFieldConfig : MagneticFieldConfig {
    std::vector<Coil> coils;

    explicit MagneticCoilsFieldConfig(std::vector<Coil> values) : coils(std::move(values)) {
    }

    void apply(Field3d& fieldB, const Domain& domain) const override;
};

// Bphi через Jz
struct BphiConfig final : MagneticFieldConfig {
    double Jz, radius, startz, endz;

    BphiConfig(double j, double r, double s, double e) : Jz(j), radius(r), startz(s), endz(e) {
    }

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

void set_Bphi(Field3d& fieldB, double Jz, double radius, double startz, double endz, const Domain& domain);
