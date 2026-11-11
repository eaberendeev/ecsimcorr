#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "grid/domain.h"
#include "utils/containers.h"

// Предписанное значение электрического DOF для Дирихле-строки.
struct ElectricDirichletEntry {
    int row;
    double value;
};

// Replace constrained equations by E_b = prescribed. Columns stay intact so
// free equations retain A_fb*E_b; RHS needs no second elimination of that term.
inline void impose_dirichlet_rows(Operator& mat, const std::vector<ElectricDirichletEntry>& entries) {
    RECORD_TIMER;
    if (entries.empty())
        return;
    mat.makeCompressed();
    double* values = mat.valuePtr();
    const int* outer = mat.outerIndexPtr();
    const int* inner = mat.innerIndexPtr();
    // A row may be prescribed by several compatible boundary conditions.
    // Give each distinct row to one thread, keeping free columns intact.
    std::vector<int> rows;
    rows.reserve(entries.size());
    for (const ElectricDirichletEntry& entry : entries) rows.push_back(entry.row);
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    int missingDiagonal = mat.rows();
#pragma omp parallel for reduction(min : missingDiagonal)
    for (size_t i = 0; i < rows.size(); ++i) {
        const int row = rows[i];
        bool hasDiag = false;
        for (int j = outer[row]; j < outer[row + 1]; ++j) {
            const bool diagonal = inner[j] == row;
            values[j] = diagonal ? 1.0 : 0.0;
            hasDiag |= diagonal;
        }
        if (!hasDiag)
            missingDiagonal = std::min(missingDiagonal, row);
    }
    if (missingDiagonal != mat.rows())
        throw std::runtime_error("Dirichlet row has no diagonal entry: " + std::to_string(missingDiagonal));
}

inline void impose_dirichlet_rhs(Field3d& rhs, const std::vector<ElectricDirichletEntry>& entries,
                                 const Domain& domain) {
    RECORD_TIMER;
    for (const auto& e : entries) {
        const auto [i, j, k, d] = domain.grid.pos_vind_range(e.row);
        rhs(i, j, k, d) = e.value;
    }
}

class FieldBoundaryCondition {
   public:
    virtual ~FieldBoundaryCondition() = default;

    FieldBoundaryCondition(Face face) : face_(face) {
    }

    // Применение к полям (например, задать значение на границе)
    virtual void apply_to_fields(Field3d& /*fields*/, FieldType /*field*/, const Domain& /*domain*/) {
        // по умолчанию ничего
    }

    // Разовая инициализация электрического поля на границе (вызывается при старте)
    virtual void init_electric_field(Field3d& /*fieldE*/, const Domain& /*domain*/) const {
        // по умолчанию ничего
    }

    virtual void apply_to_operator(Operator& /*mat*/, const Domain& /*domain*/) {
        // по умолчанию ничего
    }

    // Сбор электрических Дирихле-DOF (используется в режиме "solv" er0).
    // По умолчанию ничего: реализуют только условия, предписывающие значения поля.
    virtual void collect_electric_dirichlet(std::vector<ElectricDirichletEntry>& /*entries*/,
                                            const Domain& /*domain*/) const {
    }

    // Грань, для которой задано условие (у OpenBoundaryConditionArray — UNAVAILABLE,
    // т.к. оно может покрывать сразу несколько граней).
    virtual Face face() const {
        return face_;
    }

    Face face_;
};

// Mesh/field constraints, independent of particle transport and emission.
// Initialized before stepping; caches are rebuilt with the boundary configuration.
class FieldBoundaryConditions {
   public:
    FieldBoundaryConditions() = default;

    FieldBoundaryConditions(std::vector<std::unique_ptr<FieldBoundaryCondition>> conditions,
                            std::array<bool, 3> periodic)
        : conditions_(std::move(conditions)), periodic_(periodic) {
    }

    void apply_to_fields(Field3d& fields, FieldType field_t, const Domain& domain) const {
        RECORD_TIMER;
        for (const std::unique_ptr<FieldBoundaryCondition>& cond : conditions_) {
            cond->apply_to_fields(fields, field_t, domain);
        }
    }

    void init_electric_field(Field3d& fieldE, const Domain& domain) const {
        RECORD_TIMER;
        for (const std::unique_ptr<FieldBoundaryCondition>& cond : conditions_) {
            cond->init_electric_field(fieldE, domain);
        }
    }

    void apply_to_operator(Operator& mat, const Domain& domain) const {
        RECORD_TIMER;
        for (const std::unique_ptr<FieldBoundaryCondition>& cond : conditions_) {
            cond->apply_to_operator(mat, domain);
        }
    }

    bool is_periodic(int axis) const {
        return periodic_[axis];
    }

    // One physical periodic copy: raw [g,g+N), excluding halos and the
    // duplicate endpoint at logical N. Nonperiodic axes retain stored nodes.
    IndexRange active_range(const Grid& grid) const {
        IndexRange range;
        range.start = Vector3I(0, 0, 0);
        range.end = Vector3I(grid.size().x(), grid.size().y(), grid.size().z());

        if (is_periodic(0)) {
            range.start.x() = GHOST_CELLS;
            range.end.x() -= (GHOST_CELLS + 1);
        }
        if (is_periodic(1)) {
            range.start.y() = GHOST_CELLS;
            range.end.y() -= (GHOST_CELLS + 1);
        }
        if (is_periodic(2)) {
            range.start.z() = GHOST_CELLS;
            range.end.z() -= (GHOST_CELLS + 1);
        }
        return range;
    }

    Vector3R wrap_periodic(const Vector3R& coord, const Domain& domain) const {
        Vector3R res = coord;
        const Vector3R& box_min = domain.geom.box_min;
        const Vector3R& box_max = domain.geom.box_max;
        for (int axis = 0; axis < 3; ++axis) {
            if (!periodic_[axis])
                continue;
            const double length = box_max[axis] - box_min[axis];
            res[axis] = box_min[axis] + std::fmod(res[axis] - box_min[axis], length);
            if (res[axis] < box_min[axis])
                res[axis] += length;
            // A coordinate just below box_min rounds up to the excluded box_max.
            if (res[axis] >= box_max[axis])
                res[axis] = box_min[axis];
        }
        return res;
    }

    int wrap_index(int idx, int axis, const Domain& domain) const {
        if (!periodic_[axis])
            return idx;

        const Vector3I& size = domain.grid.size();
        // A nearest-neighbor Yee stencil can leave storage by one node.
        // Translate that neighbor by one physical period N=size-2*g-1;
        // stored halo aliases themselves keep their raw indices.
        int all_cells = size[axis] - 1;
        if (idx < 0)
            return all_cells - 2 * GHOST_CELLS + idx;
        if (idx > size[axis] - 1)
            return idx - all_cells + 2 * GHOST_CELLS;
        return idx;
    }

    const std::vector<ElectricDirichletEntry>& electric_dirichlet(const Domain& domain) const;

   protected:
    // Own each wall law once; the particle coordinator borrows only particle interfaces.
    std::vector<std::unique_ptr<FieldBoundaryCondition>> conditions_;
    std::array<bool, 3> periodic_{};

   private:
    mutable std::vector<ElectricDirichletEntry> electricDirichlet_;
    mutable bool electricDirichletReady_ = false;
};
