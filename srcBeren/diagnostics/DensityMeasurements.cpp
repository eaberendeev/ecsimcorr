#include "DensityMeasurements.h"

#include <cmath>
#include <stdexcept>

#include "ParticlesDiagnostic.h"

Field3d measure_density(const ParticlesArray& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                        ShapeType shape) {
    Field3d density(domain.size(), 1);
    if (shape == ShapeType::Linear)
        species.density_on_grid_update<ShapeType::Linear>(density);
    else if (shape == ShapeType::Quadratic)
        species.density_on_grid_update<ShapeType::Quadratic>(density);
    else
        throw std::invalid_argument("Density measurement supports Linear and Quadratic shapes");
    boundaries.apply_to_fields(density, FieldType::DENSITY, domain);
    return density;
}

Field3d measure_charge_density(const Species& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                               ShapeType shape) {
    Field3d charge(domain.size(), 1);
    charge.setZero();
    for (const auto& [name, particles] : species)
        if (!particles->is_neutral())
            charge.data() += measure_density(*particles, domain, boundaries, shape).data();
    return charge;
}

ChargeContinuityResidual measure_charge_continuity(const Field3d& before, const Field3d& after, const Field3d& current,
                                                   const Operator& divergence, double dt) {
    // Both endpoints belong to this one step, after injection and before
    // secondary emission is flushed. J is its full-step trajectory current.
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("Charge continuity requires a positive finite dt");
    const Eigen::VectorXd residual = (after.data() - before.data()) / dt + divergence * current.data();
    if (!residual.allFinite())
        throw std::runtime_error("Nonfinite charge continuity residual");
    const ChargeContinuityResidual norms{residual.norm(), residual.cwiseAbs().maxCoeff()};
    if (!std::isfinite(norms.l2) || !std::isfinite(norms.linf))
        throw std::runtime_error("Nonfinite charge continuity norm");
    return norms;
}
