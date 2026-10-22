#include <cmath>
#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "timer.h"

void Er0Condition::collect_entries(std::vector<ElectricDirichletEntry>& out, const Domain& domain) const {
    RECORD_TIMER;

    if (entriesReady_) {
        out.insert(out.end(), entries_.begin(), entries_.end());
        return;
    }

    std::vector<ElectricDirichletEntry> built;

    const auto size = domain.size();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double center_x = 0.5 * domain.num_cells().x() * dx;
    const double center_y = 0.5 * domain.num_cells().y() * dy;

    // Границы кольца зависят от выбранного профиля.
    double inner_radius = 0.0;
    double outer_radius = 0.0;
    switch (profile_) {
        case Er0Profile::Step:
            inner_radius = inner_radius_;
            outer_radius = inner_radius_ + width_;
            break;
        case Er0Profile::Cosine:
            inner_radius = r0_ - 0.5 * width_;
            outer_radius = r0_ + 0.5 * width_;
            break;
    }

    const double eps = -gap_ + eps_;
    // Кольцо задаётся только на самом граничном слое (k=ghost для ZMIN /
    // k=size_z-ghost-1 для ZMAX), остальные DOF открытого слоя получают 0.
    const int ghost = domain.grid.ghost_cells();
    const int k_face = (face_ == Face::ZMIN) ? ghost : size.z() - ghost - 1;

    constexpr double kPi = 3.14159265358979323846;

    for (int i = 0; i < size.x(); ++i) {
        for (int j = 0; j < size.y(); ++j) {
            for (int k = 0; k < size.z(); ++k) {
                for (int d = 0; d < 3; ++d) {
                    const Vector3R pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
                    // Только DOF открытого слоя этой грани.
                    if (!domain.geom.is_outside_face(face_, pos, eps))
                        continue;
                    // Уважаем остальные грани (цилиндр, боковые стенки).
                    if (!domain.geom.contains_ignoring_face(face_, pos, 0.0))
                        continue;

                    double value = 0.0;
                    if (d < 2 && k == k_face && width_ > 0.0) {
                        const double xx = pos.x() - center_x;
                        const double yy = pos.y() - center_y;
                        const double rr = std::hypot(xx, yy);
                        if (rr >= inner_radius && rr <= outer_radius && rr > 0.0) {
                            // Модуль радиального поля |Er(r)|.
                            double magnitude = 0.0;
                            if (profile_ == Er0Profile::Step) {
                                magnitude = potential_drop_ / width_;
                            } else {
                                magnitude = potential_drop_ / width_ *
                                            (1.0 + std::cos(2.0 * kPi * (rr - r0_) / width_));
                            }
                            const double coord = (d == 0) ? xx : yy;
                            value = -magnitude * coord / rr;
                        }
                    }
                    built.push_back({domain.vind(i, j, k, d), value});
                }
            }
        }
    }

    entries_ = built;
    entriesReady_ = true;
    out.insert(out.end(), entries_.begin(), entries_.end());
}

void Er0Condition::init_electric_field(Field3d& fieldE, const Domain& domain) const {
    RECORD_TIMER;

    if (face_ != Face::ZMIN && face_ != Face::ZMAX)
        return;
    if (width_ <= 0.0)
        return;

    std::vector<ElectricDirichletEntry> entries;
    collect_entries(entries, domain);

    // Пишем только ненулевые (кольцевые) значения; слой не зануляем.
    for (const auto& e : entries) {
        if (e.value == 0.0)
            continue;
        const auto [i, j, k, d] = domain.grid.pos_vind_range(e.row);
        fieldE(i, j, k, d) = e.value;
    }
}

void Er0Condition::apply_to_operator(Operator& mat, const Domain& domain) {
    if (!mode_solv_) {
        OpenBoundaryCondition::apply_to_operator(mat, domain);
        return;
    }

    RECORD_TIMER;

    mat.makeCompressed();
    double* values = mat.valuePtr();
    const int* outer = mat.outerIndexPtr();
    const int rows = 3 * domain.total_size();

    // В solv-режиме зануляем ТОЛЬКО строки открытого слоя er0-грани.
    // Столбцы не трогаем: свободные уравнения должны сохранять отклик тока
    // на предписанное граничное поле.
#pragma omp parallel for schedule(dynamic, 1024)
    for (int row = 0; row < rows; ++row) {
        const auto [i, j, k, d] = domain.grid.pos_vind_range(row);
        const Vector3R pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
        if (domain.geom.is_outside_face(face_, pos, -gap_ + eps_)) {
            for (int p = outer[row]; p < outer[row + 1]; ++p)
                values[p] = 0.0;
        }
    }
}
