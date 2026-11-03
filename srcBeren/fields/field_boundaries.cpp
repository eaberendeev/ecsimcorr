#include "field_boundaries.h"

#include "timer.h"

const std::vector<ElectricDirichletEntry>& FieldBoundaryConditions::electric_dirichlet(const Domain& domain) const {
    RECORD_TIMER;

    if (electricDirichletReady_)
        return electricDirichlet_;

    electricDirichlet_.clear();
    for (const std::unique_ptr<FieldBoundaryCondition>& cond : conditions_) {
        cond->collect_electric_dirichlet(electricDirichlet_, domain);
    }

    std::stable_sort(electricDirichlet_.begin(), electricDirichlet_.end(),
                     [](const ElectricDirichletEntry& a, const ElectricDirichletEntry& b) { return a.row < b.row; });

    for (size_t i = 1; i < electricDirichlet_.size(); ++i) {
        const auto& prev = electricDirichlet_[i - 1];
        const auto& cur = electricDirichlet_[i];
        if (prev.row == cur.row && prev.value != cur.value) {
            throw std::runtime_error("conflicting prescribed electric values on DOF row " + std::to_string(cur.row) +
                                     ": " + std::to_string(prev.value) + " and " + std::to_string(cur.value));
        }
    }

    electricDirichletReady_ = true;
    return electricDirichlet_;
}
