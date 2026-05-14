#include "external_fieldsB.h"

#include "Coil.h"
#include "Mesh.h"
#include "config.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

static void add_uniform_field(Field3d& field, double bx, double by, double bz) {
    auto sz = field.sizes();
    for (int i = 0; i < sz.x(); ++i) {
        for (int j = 0; j < sz.y(); ++j) {
            for (int k = 0; k < sz.z(); ++k) {
                field(i, j, k, Axis::X) += bx;
                field(i, j, k, Axis::Y) += by;
                field(i, j, k, Axis::Z) += bz;
            }
        }
    }
}

// Прототип локальной утилиты задания вихревого поля Bphi
void set_Bphi(Field3d& fieldB, double Jz, double radius, double startz, double endz, const Domain& domain) {
    const auto size_x = fieldB.sizes().x();
    const auto size_y = fieldB.sizes().y();
    const auto size_z = fieldB.sizes().z();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    const int sz = std::max(0, static_cast<int>(std::llround(startz / dz)));
    const int ez = std::min(static_cast<int>(size_z), static_cast<int>(std::llround(endz / dz)));

    const double center_x = 0.5 * (size_x - 3) * dx;
    const double center_y = 0.5 * (size_y - 3) * dy;
    const double eps_r = 1e-12;

    for (int k = sz; k < ez; ++k) {
        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                // X-face node
                double xx = i * dx - center_x - dx * GHOST_CELLS;
                double yy = (j + 0.5) * dy - center_y - dy * GHOST_CELLS;
                double rr = std::hypot(xx, yy);
                if (rr < radius) {
                    fieldB(i, j, k, Axis::X) = -0.5 * Jz * (yy);
                } else {
                    rr = std::max(rr, eps_r);
                    fieldB(i, j, k, Axis::X) = -0.5 * radius * radius * Jz * yy / (rr * rr);
                }

                // Y-face node
                yy = j * dy - center_y - dy * GHOST_CELLS;
                xx = (i + 0.5) * dx - center_x - dx * GHOST_CELLS;
                rr = std::hypot(xx, yy);
                if (rr < radius) {
                    fieldB(i, j, k, Axis::Y) = 0.5 * Jz * (xx);
                } else {
                    rr = std::max(rr, eps_r);
                    fieldB(i, j, k, Axis::Y) = 0.5 * radius * radius * Jz * xx / (rr * rr);
                }
            }
        }
    }
}

void MagneticUniformFieldConfig::apply(Field3d& fieldB, const Domain&) const {
    add_uniform_field(fieldB, bx, by, bz);
}

void MagneticCoilsFieldConfig::apply(Field3d& fieldB, const Domain& domain) const {
    // set_coils добавляет вклад в fieldB (использует +=), подадим конфиг в
    // нужном формате
    json cfg;
    cfg["Coils"] = coils_json;
    set_coils(fieldB, domain, cfg);
}

void BphiConfig::apply(Field3d& fieldB, const Domain& domain) const {
    set_Bphi(fieldB, Jz, radius, startz, endz, domain);
}

void MagneticCompositeFieldConfig::apply(Field3d& fieldB, const Domain& domain) const {
    for (auto& p : parts) p->apply(fieldB, domain);
}

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
            return std::make_unique<MagneticCoilsFieldConfig>(arr);
        }
    }

    // Поддержка вихревого поля Bphi из JSON: { "bphi": { "Jz": ..., "radius":
    // ..., "startz": ..., "endz": ... } } if (obj.contains("bphi")) {
    //     const auto& b = obj.at("bphi");
    //     const double Jz = b.at("Jz");
    //     const double radius = b.at("radius");
    //     const double startz = b.at("startz");
    //     const double endz = b.at("endz");
    //     return std::make_unique<BphiConfig>(Jz, radius, startz, endz);
    // }

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
