#include "external_fieldsB.h"

#include "Coil.h"

static void add_uniform_field(Field3d& field, double bx, double by, double bz) {
    auto sz = field.sizes();
    for (int i = 0; i < sz.x(); ++i) {
        for (int j = 0; j < sz.y(); ++j) {
            for (int k = 0; k < sz.z(); ++k) {
                field(i, j, k, static_cast<int>(Axis::X)) += bx;
                field(i, j, k, static_cast<int>(Axis::Y)) += by;
                field(i, j, k, static_cast<int>(Axis::Z)) += bz;
            }
        }
    }
}

// Прототип локальной утилиты задания вихревого поля Bphi
void set_Bphi(Field3d& fieldB, double Jz, double radius, double startz, double endz, const Domain& domain) {
    const auto size_x = fieldB.sizes().x();
    const auto size_y = fieldB.sizes().y();
    const auto size_z = fieldB.sizes().z();
    const Vector3R center = 0.5 * (domain.geom.box_min + domain.geom.box_max);
    const double center_x = center.x(), center_y = center.y();
    const double eps_r = 1e-12;

    // Compare physical nodes directly: conversion with ceil can misclassify an
    // endpoint that was itself computed as a Yee-node coordinate.
    for (int k = 0; k < size_z; ++k) {
        const double z = domain.get_node_position(0, 0, k, FieldType::MAGNETIC, 0).z();
        if (z < startz || z >= endz)
            continue;
        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                // X-face node
                const Vector3R bx_node = domain.get_node_position(i, j, k, FieldType::MAGNETIC, 0);
                double xx = bx_node.x() - center_x;
                double yy = bx_node.y() - center_y;
                double rr = std::hypot(xx, yy);
                if (rr < radius) {
                    fieldB(i, j, k, static_cast<int>(Axis::X)) = -0.5 * Jz * (yy);
                } else {
                    rr = std::max(rr, eps_r);
                    fieldB(i, j, k, static_cast<int>(Axis::X)) = -0.5 * radius * radius * Jz * yy / (rr * rr);
                }

                // Y-face node
                const Vector3R by_node = domain.get_node_position(i, j, k, FieldType::MAGNETIC, 1);
                yy = by_node.y() - center_y;
                xx = by_node.x() - center_x;
                rr = std::hypot(xx, yy);
                if (rr < radius) {
                    fieldB(i, j, k, static_cast<int>(Axis::Y)) = 0.5 * Jz * (xx);
                } else {
                    rr = std::max(rr, eps_r);
                    fieldB(i, j, k, static_cast<int>(Axis::Y)) = 0.5 * radius * radius * Jz * xx / (rr * rr);
                }
            }
        }
    }
}

void MagneticUniformFieldConfig::apply(Field3d& fieldB, const Domain&) const {
    add_uniform_field(fieldB, bx, by, bz);
}

void MagneticCoilsFieldConfig::apply(Field3d& fieldB, const Domain& domain) const {
    set_coils(fieldB, domain, coils);
}

void BphiConfig::apply(Field3d& fieldB, const Domain& domain) const {
    set_Bphi(fieldB, Jz, radius, startz, endz, domain);
}

void MagneticCompositeFieldConfig::apply(Field3d& fieldB, const Domain& domain) const {
    for (auto& p : parts) p->apply(fieldB, domain);
}
