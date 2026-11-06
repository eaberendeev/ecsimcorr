#include "external_fieldsE.h"

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
    const Vector3R center = 0.5 * (domain.geom.box_min + domain.geom.box_max);
    const double center_x = center.x(), center_y = center.y();
    for (auto k = 0; k < size_z; k++) {
        for (auto i = 0; i < size_x; i++) {
            for (auto j = 0; j < size_y; j++) {
                const Vector3R ex_node = domain.storage_node_position(i, j, k, FieldType::ELECTRIC, 0);
                double yy = ex_node.y() - center_y;
                double xx = ex_node.x() - center_x;
                double rr = std::hypot(xx, yy);
                // Er = 0.5*r; r < radius
                // Er = 0.5*radius*radius/r; r > radius
                if (rr == 0 || rr < r_cyl) {
                    fieldE(i, j, k, 0) = 0.5 * value * xx;   // rr * (xx / rr);
                } else {
                    fieldE(i, j, k, 0) = 0.5 * value * r_cyl * r_cyl / rr * (xx / rr);
                }

                const Vector3R ey_node = domain.storage_node_position(i, j, k, FieldType::ELECTRIC, 1);
                xx = ey_node.x() - center_x;
                yy = ey_node.y() - center_y;
                rr = std::hypot(xx, yy);
                if (rr == 0 || rr < r_cyl) {
                    fieldE(i, j, k, 1) = 0.5 * value * yy;   // rr * (yy / rr);
                } else {
                    fieldE(i, j, k, 1) = 0.5 * value * r_cyl * r_cyl / rr * (yy / rr);
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
