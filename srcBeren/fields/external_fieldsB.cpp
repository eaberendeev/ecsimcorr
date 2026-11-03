#include "external_fieldsB.h"

#include "Coil.h"
#include "config.h"

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
                    fieldB(i, j, k, static_cast<int>(Axis::X)) = -0.5 * Jz * (yy);
                } else {
                    rr = std::max(rr, eps_r);
                    fieldB(i, j, k, static_cast<int>(Axis::X)) = -0.5 * radius * radius * Jz * yy / (rr * rr);
                }

                // Y-face node
                yy = j * dy - center_y - dy * GHOST_CELLS;
                xx = (i + 0.5) * dx - center_x - dx * GHOST_CELLS;
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
