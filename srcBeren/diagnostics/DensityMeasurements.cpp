#include "DensityMeasurements.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "coupling/density.h"
#include "runtime/profiler.h"
#include "utils/sparse.h"
#include "utils/thread_partitioned_matrix.h"

// Unfolded deposit: the caller applies the density boundary condition.
static void deposit_density(const ParticlesArray& species, ShapeType shape, Field3d& density) {
    if (shape == ShapeType::Linear)
        coupling::deposit_density<ShapeType::Linear>(species, density);
    else if (shape == ShapeType::Quadratic)
        coupling::deposit_density<ShapeType::Quadratic>(species, density);
    else
        throw std::invalid_argument("Density measurement supports Linear and Quadratic shapes");
}

Field3d measure_density(const ParticlesArray& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                        ShapeType shape) {
    RECORD_TIMER;
    Field3d density(domain.size(), 1);
    deposit_density(species, shape, density);
    boundaries.apply_to_fields(density, FieldType::DENSITY, domain);
    return density;
}

Field3d measure_charge_density(const Species& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                               ShapeType shape) {
    RECORD_TIMER;
    Field3d charge = Field3d::Zero(domain.size(), 1);
    Field3d density(domain.size(), 1);
    for (const Species::value_type& entry : species) {
        if (entry.second->is_neutral())
            continue;
        deposit_density(*entry.second, shape, density);
        charge += density;
    }
    // The density conditions are linear: one application to the sum serves every species.
    boundaries.apply_to_fields(charge, FieldType::DENSITY, domain);
    return charge;
}

template <typename Divergence>
static ChargeContinuityResidual measure_charge_continuity_impl(const Field3d& before, const Field3d& after,
                                                               const Field3d& current, const Divergence& divergence,
                                                               double dt) {
    RECORD_TIMER;
    // Both endpoints belong to this one step, after injection and before
    // secondary emission is flushed. J is its full-step trajectory current.
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("Charge continuity requires a positive finite dt");
    Field3d residual = (after - before) * (1.0 / dt);
    // div maps a three-component current to a scalar field. Supply the scalar
    // result explicitly: matrix * Field3d otherwise keeps the input dimensions.
    Field3d div_current(before.sizes(), 1);
    spmv(divergence, current, div_current);
    residual += div_current;
    bool nonfinite = false;
    double linf = 0;
#pragma omp parallel for reduction(|| : nonfinite) reduction(max : linf)
    for (int i = 0; i < residual.capacity(); ++i) {
        const double value = residual[i];
        nonfinite = nonfinite || !std::isfinite(value);
        linf = std::max(linf, std::abs(value));
    }
    if (nonfinite)
        throw std::runtime_error("Nonfinite charge continuity residual");
    const ChargeContinuityResidual norms{residual.norm(), linf};
    if (!std::isfinite(norms.l2) || !std::isfinite(norms.linf))
        throw std::runtime_error("Nonfinite charge continuity norm");
    return norms;
}

ChargeContinuityResidual measure_charge_continuity(const Field3d& before, const Field3d& after, const Field3d& current,
                                                   const Operator& divergence, double dt) {
    return measure_charge_continuity_impl(before, after, current, divergence, dt);
}

ChargeContinuityResidual measure_charge_continuity(const Field3d& before, const Field3d& after, const Field3d& current,
                                                   const ThreadPartitionedSparseMatrix<double>& divergence, double dt) {
    return measure_charge_continuity_impl(before, after, current, divergence, dt);
}
