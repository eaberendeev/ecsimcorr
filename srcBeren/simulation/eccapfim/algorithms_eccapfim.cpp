#include "algorithms_eccapfim.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "Shape.h"
#include "decompose_esirkepov_current.h"

namespace algorithmsECCAPFIM {
int periodic_node(int index, int axis, const Domain& domain) {
    return domain.grid.periodic_node_owner(index + domain.grid.ghost_cells(), axis);
}

Vector3R wrap_position(const Vector3R& position, const Domain& domain) {
    Vector3R result = position;
    for (int d = 0; d < 3; ++d) {
        const double length = domain.num_cells(d) * domain.cell_size(d);
        result[d] -= length * std::floor(result[d] / length);
        if (result[d] >= length)
            result[d] = 0;
    }
    return result;
}

void copy_periodic_nodes(Field3d& field, const Domain& domain) {
    const auto size = domain.size();
    const int ghost = domain.grid.ghost_cells();
    // Copy from the unique N^3 physical nodes; never add already folded sources.
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int k = 0; k < size.z(); ++k)
                for (int c = 0; c < field.nd(); ++c)
                    field(i, j, k, c) = field(periodic_node(i - ghost, 0, domain), periodic_node(j - ghost, 1, domain),
                                              periodic_node(k - ghost, 2, domain), c);
}

namespace {
std::vector<double> orbit_cuts(const Vector3R& start, const Vector3R& end) {
    std::vector<double> cuts{0, 1};
    for (int d = 0; d < 3; ++d) {
        if (!std::isfinite(start[d]) || !std::isfinite(end[d]))
            throw std::runtime_error("ECCAPFIM: nonfinite particle orbit");
        const double lo = std::min(start[d], end[d]), hi = std::max(start[d], end[d]);
        if (std::abs(lo) > 1e8 || std::abs(hi) > 1e8 || hi - lo > 100000)
            throw std::runtime_error("ECCAPFIM: orbit exceeds traversal limit");
        for (double plane = std::floor(lo - 0.5) + 1.5; plane < hi; plane += 1) {
            const double t = (plane - start[d]) / (end[d] - start[d]);
            if (t > 0 && t < 1)
                cuts.push_back(t);
        }
    }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
}

Vector3I fill_orbit_weights(const Vector3R& start, const Vector3R& end, CurrentBuffer<3>& weights) {
    const Vector3R midpoint = 0.5 * (start + end);
    Vector3I center;
    for (int d = 0; d < 3; ++d) center[d] = static_cast<int>(std::floor(midpoint[d] + 0.5));

    ParticleShape<Shape2, 3> start_shape, end_shape;
    ParticleShape<Shape, 2> midpoint_shape;
    // Both endpoints belong to the same quadratic spline interval after
    // splitting at half-cell planes. All shapes share start_ = center - 1.
    start_shape.fill_from_normalized(start, center, 1);
    end_shape.fill_from_normalized(end, center, 1);
    midpoint_shape.fill_from_normalized(midpoint - Vector3R(0.5), center, 1);
    weights.zero();
    for (int n = 0; n < 3; ++n)
        for (int m = 0; m < 3; ++m)
            for (int k = 0; k < 3; ++k) {
                // Esirkepov transverse factors with a linear longitudinal
                // factor at the midpoint: quadratic charge -> edge current.
                if (n < 2)
                    weights(n, m, k, X) = midpoint_shape(n, X) *
                                          (end_shape(m, Y) * (2 * end_shape(k, Z) + start_shape(k, Z)) +
                                           start_shape(m, Y) * (2 * start_shape(k, Z) + end_shape(k, Z))) /
                                          6;
                if (m < 2)
                    weights(n, m, k, Y) = midpoint_shape(m, Y) *
                                          (end_shape(n, X) * (2 * end_shape(k, Z) + start_shape(k, Z)) +
                                           start_shape(n, X) * (2 * start_shape(k, Z) + end_shape(k, Z))) /
                                          6;
                if (k < 2)
                    weights(n, m, k, Z) = midpoint_shape(k, Z) *
                                          (end_shape(n, X) * (2 * end_shape(m, Y) + start_shape(m, Y)) +
                                           start_shape(n, X) * (2 * start_shape(m, Y) + end_shape(m, Y))) /
                                          6;
            }
    return start_shape.start_;
}

Vector3R gather_magnetic(const Field3d& field, const Vector3R& coord, const Domain& domain) {
    const auto r = domain.to_cell_coordinates(coord);
    Vector3R result{0, 0, 0};
    for (int c = 0; c < 3; ++c) {
        const auto shift = magnetic_shift(c);
        const auto shifted = r - shift;
        Vector3I center;
        for (int d = 0; d < 3; ++d) center[d] = static_cast<int>(std::floor(shifted[d] + 0.5));
        ParticleShape<Shape2, 3> shape;
        shape.fill_from_normalized(shifted, center, 1);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) {
                    const auto index = shape.start_ + Vector3I(i, j, k);
                    const double w = shape(i, X) * shape(j, Y) * shape(k, Z);
                    result[c] += w * field(periodic_node(index.x(), 0, domain), periodic_node(index.y(), 1, domain),
                                           periodic_node(index.z(), 2, domain), c);
                }
    }
    return result;
}
}   // namespace

OrbitFields gather_orbit(const Field3d& electric, const Field3d& magnetic, const Vector3R& start, const Vector3R& end,
                         const Domain& domain) {
    OrbitFields result;
    const auto a = domain.to_cell_coordinates(start), b = domain.to_cell_coordinates(end);
    const auto cuts = orbit_cuts(a, b);
    CurrentBuffer<3> weights;
    for (size_t s = 1; s < cuts.size(); ++s) {
        const Vector3R r0 = a + cuts[s - 1] * (b - a), rn = a + cuts[s] * (b - a);
        const Vector3I first = fill_orbit_weights(r0, rn, weights);
        const double fraction = cuts[s] - cuts[s - 1];
        for (int n = 0; n < 3; ++n) {
            const int ix = periodic_node(first.x() + n, X, domain);
            for (int m = 0; m < 3; ++m) {
                const int iy = periodic_node(first.y() + m, Y, domain);
                for (int k = 0; k < 3; ++k) {
                    const int iz = periodic_node(first.z() + k, Z, domain);
                    for (int c = 0; c < 3; ++c)
                        result.electric[c] += fraction * weights(n, m, k, c) * electric(ix, iy, iz, c);
                }
            }
        }
        const double t = 0.5 * (cuts[s] + cuts[s - 1]);
        result.magnetic += fraction * gather_magnetic(magnetic, start + t * (end - start), domain);
    }
    return result;
}

void deposit_orbit(Field3d& current, const Vector3R& start, const Vector3R& end, const Vector3R& velocity,
                   double charge_weight, const Domain& domain) {
    const auto a = domain.to_cell_coordinates(start), b = domain.to_cell_coordinates(end);
    const auto cuts = orbit_cuts(a, b);
    CurrentBuffer<3> weights;
    for (size_t s = 1; s < cuts.size(); ++s) {
        const Vector3R r0 = a + cuts[s - 1] * (b - a), rn = a + cuts[s] * (b - a);
        const Vector3I first = fill_orbit_weights(r0, rn, weights);
        const Vector3R value = (charge_weight * (cuts[s] - cuts[s - 1])) * velocity;
        for (int n = 0; n < 3; ++n) {
            const int ix = periodic_node(first.x() + n, X, domain);
            for (int m = 0; m < 3; ++m) {
                const int iy = periodic_node(first.y() + m, Y, domain);
                for (int k = 0; k < 3; ++k) {
                    const int iz = periodic_node(first.z() + k, Z, domain);
                    for (int c = 0; c < 3; ++c) {
#pragma omp atomic update
                        current(ix, iy, iz, c) += value[c] * weights(n, m, k, c);
                    }
                }
            }
        }
    }
}

void deposit_charge(Field3d& density, const Vector3R& position, double charge_weight, const Domain& domain) {
    const auto r = domain.to_cell_coordinates(wrap_position(position, domain));
    Vector3I center;
    for (int d = 0; d < 3; ++d) center[d] = static_cast<int>(std::floor(r[d] + 0.5));
    ParticleShape<Shape2, 3> shape;
    shape.fill_from_normalized(r, center, 1);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                const auto index = shape.start_ + Vector3I(i, j, k);
                const double value = charge_weight * shape(i, X) * shape(j, Y) * shape(k, Z);
#pragma omp atomic update
                density(periodic_node(index.x(), 0, domain), periodic_node(index.y(), 1, domain),
                        periodic_node(index.z(), 2, domain), 0) += value;
            }
}

Particle advance_particle(const Particle& initial, const Field3d& electric, const Field3d& magnetic, Field3d& current,
                          double charge_mass, double charge_weight, double dt, const Domain& domain, double tolerance,
                          int max_iterations) {
    if (!(dt > 0) || !std::isfinite(dt) || !(tolerance > 0) || !std::isfinite(tolerance) || max_iterations < 1)
        throw std::invalid_argument("ECCAPFIM: invalid particle solver parameters");
    for (int d = 0; d < 3; ++d)
        if (!std::isfinite(initial.coord[d]) || !std::isfinite(initial.velocity[d]))
            throw std::runtime_error("ECCAPFIM: nonfinite initial particle");
    if (charge_mass == 0) {
        Particle result = initial;
        result.move(dt);
        for (int d = 0; d < 3; ++d)
            if (!std::isfinite(result.coord[d]))
                throw std::runtime_error("ECCAPFIM: nonfinite neutral orbit");
        result.coord = wrap_position(result.coord, domain);
        return result;
    }
    Particle start = initial, end = initial;
    double remaining = dt;
    int substeps = 0;
    while (remaining > 0) {
        if (++substeps > 10000)
            throw std::runtime_error("ECCAPFIM: too many particle substeps");
        double tau = remaining;
        bool converged = false;
        for (int retry = 0; retry < 30 && !converged; ++retry) {
            end = start;
            end.coord += tau * start.velocity;
            for (int it = 0; it < max_iterations; ++it) {
                const auto fields = gather_orbit(electric, magnetic, start.coord, end.coord, domain);
                const auto b = (0.5 * tau * charge_mass) * fields.magnetic;
                const auto w = start.velocity + (0.5 * tau * charge_mass) * fields.electric;
                const auto vh = (w + w.cross(b) + b * w.dot(b)) / (1 + b.dot(b));
                const auto next_position = start.coord + tau * vh;
                const auto next_velocity = 2.0 * vh - start.velocity;
                const double error =
                    (next_position - end.coord).length() + tau * (next_velocity - end.velocity).length();
                end.coord = next_position;
                end.velocity = next_velocity;
                if (!std::isfinite(error))
                    throw std::runtime_error("ECCAPFIM: nonfinite particle iterate");
                if (error <= tolerance * (1 + tau * vh.length())) {
                    // Check the Lorentz equation using fields on the accepted orbit.
                    const auto accepted = gather_orbit(electric, magnetic, start.coord, end.coord, domain);
                    const auto residual = end.velocity - start.velocity -
                                          (tau * charge_mass) * (accepted.electric + vh.cross(accepted.magnetic));
                    converged = residual.length() <= tolerance * (1 + end.velocity.length());
                    if (converged)
                        break;
                }
            }
            if (!converged)
                tau *= 0.5;
        }
        if (!converged || !(tau > 0) || remaining - tau == remaining)
            throw std::runtime_error("ECCAPFIM: particle solve did not converge");
        const auto vh = 0.5 * (start.velocity + end.velocity);
        deposit_orbit(current, start.coord, end.coord, vh, charge_weight * (tau / dt), domain);
        remaining = tau >= remaining ? 0 : remaining - tau;
        start = end;
        // Keep coordinates bounded after depositing each unwrapped substep.
        start.coord = wrap_position(start.coord, domain);
    }
    return start;
}
}   // namespace algorithmsECCAPFIM
