#include "algorithms_eccapfim.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "coupling/shape.h"
#include "coupling/esirkepov.h"

namespace {
// Historical axis spellings used by the XPIC kernels.
constexpr int X = static_cast<int>(Axis::X);
constexpr int Y = static_cast<int>(Axis::Y);
constexpr int Z = static_cast<int>(Axis::Z);
}   // namespace

namespace algorithmsECCAPFIM {
int periodic_node(int index, int axis, const Domain& domain) {
    return domain.grid.periodic_node_owner(index + domain.grid.ghost_cells(), axis);
}

int node_owner(int index, int axis, const Domain& domain, const BoundaryConditionHandler* bc) {
    if (!bc || bc->is_periodic(axis))
        return periodic_node(index, axis, domain);
    const int stored = index + domain.grid.ghost_cells();
    if (stored < 0 || stored >= domain.size()[axis])
        throw std::runtime_error("ECCAPFIM: orbit stencil outside nonperiodic storage");
    return stored;
}

Vector3R wrap_position(const Vector3R& position, const Domain& domain, const BoundaryConditionHandler* bc) {
    if (bc)
        return bc->wrap_periodic(position, domain);
    Vector3R result = position;
    for (int d = 0; d < 3; ++d) {
        const double length = domain.num_cells(d) * domain.cell_size(d);
        result[d] -= length * std::floor(result[d] / length);
        if (result[d] >= length)
            result[d] = 0;
    }
    return result;
}

void copy_periodic_nodes(Field3d& field, const Domain& domain, const BoundaryConditionHandler* bc) {
    const auto size = domain.size();
    const int ghost = domain.grid.ghost_cells();
    // Copy from the unique N^3 physical nodes; never add already folded sources.
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int k = 0; k < size.z(); ++k)
                for (int c = 0; c < field.nd(); ++c)
                    field(i, j, k, c) = field(node_owner(i - ghost, 0, domain, bc), node_owner(j - ghost, 1, domain, bc),
                                              node_owner(k - ghost, 2, domain, bc), c);
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

Vector3R gather_magnetic(const Field3d& field, const Vector3R& coord, const Domain& domain, const BoundaryConditionHandler* bc) {
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
                    result[c] += w * field(node_owner(index.x(), 0, domain, bc), node_owner(index.y(), 1, domain, bc),
                                           node_owner(index.z(), 2, domain, bc), c);
                }
    }
    return result;
}
}   // namespace

OrbitFields gather_orbit(const Field3d& electric, const Field3d& magnetic, const Vector3R& start, const Vector3R& end,
                         const Domain& domain, const BoundaryConditionHandler* bc) {
    OrbitFields result;
    const auto a = domain.to_cell_coordinates(start), b = domain.to_cell_coordinates(end);
    const auto cuts = orbit_cuts(a, b);
    CurrentBuffer<3> weights;
    for (size_t s = 1; s < cuts.size(); ++s) {
        const Vector3R r0 = a + cuts[s - 1] * (b - a), rn = a + cuts[s] * (b - a);
        const Vector3I first = fill_orbit_weights(r0, rn, weights);
        const double fraction = cuts[s] - cuts[s - 1];
        for (int n = 0; n < 3; ++n) {
            const int ix = node_owner(first.x() + n, X, domain, bc);
            for (int m = 0; m < 3; ++m) {
                const int iy = node_owner(first.y() + m, Y, domain, bc);
                for (int k = 0; k < 3; ++k) {
                    const int iz = node_owner(first.z() + k, Z, domain, bc);
                    for (int c = 0; c < 3; ++c)
                        result.electric[c] += fraction * weights(n, m, k, c) * electric(ix, iy, iz, c);
                }
            }
        }
        const double t = 0.5 * (cuts[s] + cuts[s - 1]);
        result.magnetic += fraction * gather_magnetic(magnetic, start + t * (end - start), domain, bc);
    }
    return result;
}

void deposit_orbit(Field3d& current, const Vector3R& start, const Vector3R& end, const Vector3R& velocity,
                   double charge_weight, const Domain& domain, const BoundaryConditionHandler* bc) {
    const auto a = domain.to_cell_coordinates(start), b = domain.to_cell_coordinates(end);
    const auto cuts = orbit_cuts(a, b);
    CurrentBuffer<3> weights;
    for (size_t s = 1; s < cuts.size(); ++s) {
        const Vector3R r0 = a + cuts[s - 1] * (b - a), rn = a + cuts[s] * (b - a);
        const Vector3I first = fill_orbit_weights(r0, rn, weights);
        const Vector3R value = (charge_weight * (cuts[s] - cuts[s - 1])) * velocity;
        for (int n = 0; n < 3; ++n) {
            const int ix = node_owner(first.x() + n, X, domain, bc);
            for (int m = 0; m < 3; ++m) {
                const int iy = node_owner(first.y() + m, Y, domain, bc);
                for (int k = 0; k < 3; ++k) {
                    const int iz = node_owner(first.z() + k, Z, domain, bc);
                    for (int c = 0; c < 3; ++c) {
#pragma omp atomic update
                        current(ix, iy, iz, c) += value[c] * weights(n, m, k, c);
                    }
                }
            }
        }
    }
}

void deposit_charge(Field3d& density, const Vector3R& position, double charge_weight, const Domain& domain, const BoundaryConditionHandler* bc) {
    const auto r = domain.to_cell_coordinates(wrap_position(position, domain, bc));
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
                density(node_owner(index.x(), 0, domain, bc), node_owner(index.y(), 1, domain, bc),
                        node_owner(index.z(), 2, domain, bc), 0) += value;
            }
}

Particle advance_particle(const Particle& initial, const Field3d& electric, const Field3d& magnetic, Field3d& current,
                          double charge_mass, double charge_weight, double dt, const Domain& domain, double tolerance,
                          int max_iterations, const BoundaryConditionHandler* bc, const PathPolicy& policy,
                          OrbitBoundaryEvents* events) {
    if (!(dt > 0) || !std::isfinite(dt) || !(tolerance > 0) || !std::isfinite(tolerance) || max_iterations < 1)
        throw std::invalid_argument("ECCAPFIM: invalid particle solver parameters");
    for (int d = 0; d < 3; ++d)
        if (!std::isfinite(initial.coord[d]) || !std::isfinite(initial.velocity[d]))
            throw std::runtime_error("ECCAPFIM: nonfinite initial particle");
    OrbitBoundaryEvents local_events;
    if (!events)
        events = &local_events;
    *events = {};
    Particle start = initial, end = initial;
    double remaining = dt;
    int substeps = 0;
    while (remaining > 0) {
        if (++substeps > 10000)
            throw std::runtime_error("ECCAPFIM: too many particle substeps");
        double tau = remaining;
        bool converged = false;
        BoundaryCrossing hit;
        // A particle exactly on the wall and moving out interacts at tau=0.
        Vector3R clipped;
        hit = shape_particle_path(start.coord, start.coord + tau * start.velocity, domain, policy, clipped);
        if (!(hit.clipped && hit.t_hit == 0)) {
            hit = {};
            for (int retry = 0; retry < 60 && !converged; ++retry) {
                end = start;
                end.coord += tau * start.velocity;
                for (int it = 0; it < max_iterations; ++it) {
                    // Trial orbits can cross a wall. Gather only up to that wall;
                    // shorten tau below and solve the Lorentz equation again.
                    shape_particle_path(start.coord, end.coord, domain, policy, clipped);
                    const auto fields = gather_orbit(electric, magnetic, start.coord, clipped, domain, bc);
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
                        converged = true;
                        break;
                    }
                }
                const auto crossing = shape_particle_path(start.coord, end.coord, domain, policy, clipped);
                if (converged && crossing.clipped) {
                    if ((end.coord - clipped).length() <= tolerance * (1 + end.coord.length())) {
                        end.coord = clipped;
                        hit = crossing;
                    } else {
                        tau *= crossing.t_hit;
                        converged = false;
                        continue;
                    }
                }
                if (converged && !hit.clipped) {
                    // Roundoff may leave a wall-limited orbit just inside the
                    // box. Snap only an outward endpoint within orbit tolerance.
                    const Face low[] = {Face::XMIN, Face::YMIN, Face::ZMIN};
                    const Face high[] = {Face::XMAX, Face::YMAX, Face::ZMAX};
                    const auto vh = 0.5 * (start.velocity + end.velocity);
                    for (int d = 0; d < 3; ++d) {
                        const Face face = vh[d] < 0 ? low[d] : high[d];
                        const double wall = vh[d] < 0 ? domain.geom.box_min[d] : domain.geom.box_max[d];
                        if (vh[d] != 0 && policy.clips(face) &&
                            std::abs(end.coord[d] - wall) <= tolerance * (1 + end.coord.length())) {
                            snap_path_end_to_plane(face, domain, end.coord);
                            hit = {true, face, 1.0, end.coord};
                            break;
                        }
                    }
                }
                if (converged) {
                    const auto accepted = gather_orbit(electric, magnetic, start.coord, end.coord, domain, bc);
                    const auto vh = 0.5 * (start.velocity + end.velocity);
                    const auto residual = end.velocity - start.velocity -
                                          (tau * charge_mass) * (accepted.electric + vh.cross(accepted.magnetic));
                    converged = residual.length() <= tolerance * (1 + end.velocity.length());
                }
                if (converged)
                    break;
                hit = {};
                tau *= 0.5;
            }
            if (!converged || !(tau > 0) || remaining - tau == remaining)
                throw std::runtime_error("ECCAPFIM: particle solve did not converge");
            const auto vh = 0.5 * (start.velocity + end.velocity);
            if (charge_weight != 0)
                deposit_orbit(current, start.coord, end.coord, vh, charge_weight * (tau / dt), domain, bc);
            remaining = tau >= remaining ? 0 : remaining - tau;
            start = end;
        }
        if (hit.clipped) {
            start.coord = wrap_position(hit.point, domain, bc);
            const auto fate = policy.classify ? policy.classify(start, hit.face)
                                              : ParticleFateResult{ParticleFate::Removed, hit.face};
            if (fate.fate == ParticleFate::Removed) {
                events->removed = hit.face;
                return start;
            }
            if (fate.fate != ParticleFate::Reflected)
                throw std::runtime_error("ECCAPFIM: wall crossing has no particle fate");
            // Specular reflection does no work. Continue the remaining physical
            // time with a new implicit orbit and the same full-dt current scale.
            const int axis = (hit.face == Face::XMIN || hit.face == Face::XMAX) ? X
                             : (hit.face == Face::YMIN || hit.face == Face::YMAX) ? Y : Z;
            start.velocity[axis] = -start.velocity[axis];
            ++events->reflected[static_cast<size_t>(hit.face)];
        }
        start.coord = wrap_position(start.coord, domain, bc);
    }
    // Box containment is half-open; an exact reflected high-wall endpoint
    // belongs to the adjacent interior cell when rebinned.
    if (bc)
        for (int d = 0; d < 3; ++d)
            if (!bc->is_periodic(d) && start.coord[d] == domain.geom.box_max[d])
                start.coord[d] = std::nextafter(start.coord[d], domain.geom.box_min[d]);
    return start;
}
}   // namespace algorithmsECCAPFIM
