#include "external_fieldsE.h"

#include "Mesh.h"
#include "config.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

// Локальные хелперы для сложения полей
static void add_uniform_field(Field3d& field, double ex, double ey, double ez) {
    auto sz = field.sizes();
    for (int i = 0; i < sz.x(); ++i) {
        for (int j = 0; j < sz.y(); ++j) {
            for (int k = 0; k < sz.z(); ++k) {
                field(i, j, k, Axis::X) += ex;
                field(i, j, k, Axis::Y) += ey;
                field(i, j, k, Axis::Z) += ez;
            }
        }
    }
}

// Сложение dst += src
static void add_field(Field3d& dst, const Field3d& src) {
    auto sz = dst.sizes();
    for (int i = 0; i < sz.x(); ++i) {
        for (int j = 0; j < sz.y(); ++j) {
            for (int k = 0; k < sz.z(); ++k) {
                dst(i, j, k, Axis::X) += src(i, j, k, Axis::X);
                dst(i, j, k, Axis::Y) += src(i, j, k, Axis::Y);
                dst(i, j, k, Axis::Z) += src(i, j, k, Axis::Z);
            }
        }
    }
}

void ElectricUniformFieldConfig::apply(Field3d& fieldE, const Domain&) const {
    // Однородный вклад добавляется, а не затирается
    add_uniform_field(fieldE, ex, ey, ez);
}

void set_uniformly_charged_cylinder(Field3d& fieldE, const Domain& domain, const double r_cyl, const double value) {
    const int size_x = fieldE.sizes().x();
    const int size_y = fieldE.sizes().y();
    const int size_z = fieldE.sizes().z();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double center_x = 0.5 * (size_x - 3) * dx;
    const double center_y = 0.5 * (size_y - 3) * dy;
    for (auto k = 0; k < size_z; k++) {
        for (auto i = 0; i < size_x; i++) {
            for (auto j = 0; j < size_y; j++) {
                double yy = j * dy - center_y - dy * GHOST_CELLS;
                double xx = (i + 0.5) * dx - center_x - dx * GHOST_CELLS;
                double rr = std::hypot(xx, yy);
                // Er = 0.5*r; r < radius
                // Er = 0.5*radius*radius/r; r > radius
                if (rr < r_cyl) {
                    fieldE(i, j, k, 0) = 0.5 * value * xx;   // rr * (xx / rr);
                } else if (rr <= 0.5 * dx * size_x) {
                    fieldE(i, j, k, 0) = 0.5 * value * r_cyl * r_cyl / rr * (xx / rr);
                } else {
                    fieldE(i, j, k, 0) = 0.;
                }

                xx = i * dx - center_x - dx * GHOST_CELLS;
                yy = (j + 0.5) * dy - center_y - dy * GHOST_CELLS;
                rr = std::hypot(xx, yy);
                if (rr < r_cyl) {
                    fieldE(i, j, k, 1) = 0.5 * value * yy;   // rr * (yy / rr);
                } else if (rr <= 0.5 * dx * size_x) {
                    fieldE(i, j, k, 1) = 0.5 * value * r_cyl * r_cyl / rr * (yy / rr);
                } else {
                    fieldE(i, j, k, 1) = 0.;
                }
            }
        }
    }
}

void ElectricUniformlyChargedCylinderConfig::apply(Field3d& fieldE, const Domain& domain) const {
    // Т.к. set_uniformly_charged_cylinder перезаписывает компоненты,
    // формируем временное поле и потом добавляем
    Field3d tmp = fieldE;   // копия по геометрии
    tmp.setZero();
    set_uniformly_charged_cylinder(tmp, domain, radius, value);
    add_field(fieldE, tmp);
}

void ElectricCompositeFieldConfig::apply(Field3d& fieldE, const Domain& domain) const {
    for (auto& p : parts) p->apply(fieldE, domain);
}

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
