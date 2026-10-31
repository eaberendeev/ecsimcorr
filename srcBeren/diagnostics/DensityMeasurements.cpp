#include "DensityMeasurements.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "coupling/density.h"
#include "sparse.h"
#include "timer.h"

Field3d measure_density(const ParticlesArray& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                        ShapeType shape) {
    RECORD_TIMER;
    Field3d density(domain.size(), 1);
    if (shape == ShapeType::Linear)
        coupling::deposit_density<ShapeType::Linear>(species, density);
    else if (shape == ShapeType::Quadratic)
        coupling::deposit_density<ShapeType::Quadratic>(species, density);
    else
        throw std::invalid_argument("Density measurement supports Linear and Quadratic shapes");
    boundaries.apply_to_fields(density, FieldType::DENSITY, domain);
    return density;
}

Field3d measure_charge_density(const Species& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                               ShapeType shape) {
    RECORD_TIMER;
    Field3d charge(domain.size(), 1);
    charge.setZero();
    for (const auto& [name, particles] : species)
        if (!particles->is_neutral())
            charge += measure_density(*particles, domain, boundaries, shape);
    return charge;
}

ChargeContinuityResidual measure_charge_continuity(const Field3d& before, const Field3d& after, const Field3d& current,
                                                   const Operator& divergence, double dt) {
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
