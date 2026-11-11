#include "algorithms_eccapfim.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "coupling/stencil.h"
#include "coupling/esirkepov.h"
#include "runtime/profiler.h"

namespace {
// Historical axis spellings used by the XPIC kernels.
constexpr int X = static_cast<int>(Axis::X);
constexpr int Y = static_cast<int>(Axis::Y);
constexpr int Z = static_cast<int>(Axis::Z);
}   // namespace

namespace algorithmsECCAPFIM {
int periodic_node(int index, int axis, const Domain& domain) {
    return domain.grid.periodic_node_owner(index + GHOST_CELLS, axis);
}

int node_owner(int index, int axis, const Domain& domain, const BoundaryConditionHandler* bc) {
    if (!bc || bc->is_periodic(axis))
        return periodic_node(index, axis, domain);
    const int stored = index + GHOST_CELLS;
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
    RECORD_TIMER;
    const auto size = domain.size();
    const int ghost = GHOST_CELLS;
    // Copy from the unique N^3 physical nodes; never add already folded sources.
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int k = 0; k < size.z(); ++k)
                for (int c = 0; c < field.nd(); ++c)
                    field(i, j, k, c) =
                        field(node_owner(i - ghost, 0, domain, bc), node_owner(j - ghost, 1, domain, bc),
                              node_owner(k - ghost, 2, domain, bc), c);
}

namespace {
// Anderson acceleration (history 2) of the fixed point v_h = G(v_h) of one
// Crank-Nicolson orbit, G(v_h) = R(B)(v^n + alpha E) with the fields gathered
// along x^n -> x^n + tau v_h. With omega_c tau >> 1 the E x B displacement
// responds to the field gradient with a factor ~ tau |grad E| / B, and plain
// Picard iteration can stall or cycle when the orbit end sits on a cell
// plane, where the orbit weights have a kink. The fixed point (the CN orbit)
// is unchanged; only the iterates are combined.
class OrbitAnderson {
   public:
    // Next iterate from the current one `x` and its image `g` = G(x).
    Vector3R next(const Vector3R& x, const Vector3R& g) {
        const Vector3R f = g - x;
        Vector3R result = g;
        if (count_ == 1) {
            const Vector3R dF = f - f_[0], dG = g - g_[0];
            const double denominator = dF.dot(dF);
            if (denominator > 1e-300)
                result = g - (dF.dot(f) / denominator) * dG;
        } else if (count_ == 2) {
            // gamma = argmin |f - dF gamma| over the two previous differences.
            const Vector3R dF0 = f - f_[0], dF1 = f - f_[1];
            const double a00 = dF0.dot(dF0), a01 = dF0.dot(dF1), a11 = dF1.dot(dF1);
            const double b0 = dF0.dot(f), b1 = dF1.dot(f);
            const double det = a00 * a11 - a01 * a01;
            if (det > 1e-12 * a00 * a11 && det > 1e-300) {
                const double gamma0 = (b0 * a11 - b1 * a01) / det, gamma1 = (b1 * a00 - b0 * a01) / det;
                result = g - gamma0 * (g - g_[0]) - gamma1 * (g - g_[1]);
            } else if (a00 > 1e-300) {
                result = g - (b0 / a00) * (g - g_[0]);
            }
        }
        if (!std::isfinite(result.x()) || !std::isfinite(result.y()) || !std::isfinite(result.z()))
            result = g;
        // Keep the newest pair first.
        f_[1] = f_[0];
        g_[1] = g_[0];
        f_[0] = f;
        g_[0] = g;
        count_ = std::min(count_ + 1, 2);
        return result;
    }

   private:
    Vector3R f_[2], g_[2];
    int count_ = 0;
};

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

    const Vector3I origin = center - Vector3I(1);
    const Stencil3D<3> start_shape = make_stencil_at<Shape2, 3>(start, origin);
    const Stencil3D<3> end_shape = make_stencil_at<Shape2, 3>(end, origin);
    const Stencil3D<2> midpoint_shape = make_stencil_at<Shape, 2>(midpoint - Vector3R(0.5), origin);
    weights.zero();
    for (int n = 0; n < 3; ++n)
        for (int m = 0; m < 3; ++m)
            for (int k = 0; k < 3; ++k) {
                // Esirkepov transverse factors with a linear longitudinal
                // factor at the midpoint: quadratic charge -> edge current.
                if (n < 2)
                    weights(n, m, k, X) = midpoint_shape(n, Axis::X) *
                                          (end_shape(m, Axis::Y) *
                                               (2 * end_shape(k, Axis::Z) + start_shape(k, Axis::Z)) +
                                           start_shape(m, Axis::Y) *
                                               (2 * start_shape(k, Axis::Z) + end_shape(k, Axis::Z))) /
                                          6;
                if (m < 2)
                    weights(n, m, k, Y) = midpoint_shape(m, Axis::Y) *
                                          (end_shape(n, Axis::X) *
                                               (2 * end_shape(k, Axis::Z) + start_shape(k, Axis::Z)) +
                                           start_shape(n, Axis::X) *
                                               (2 * start_shape(k, Axis::Z) + end_shape(k, Axis::Z))) /
                                          6;
                if (k < 2)
                    weights(n, m, k, Z) = midpoint_shape(k, Axis::Z) *
                                          (end_shape(n, Axis::X) *
                                               (2 * end_shape(m, Axis::Y) + start_shape(m, Axis::Y)) +
                                           start_shape(n, Axis::X) *
                                               (2 * start_shape(m, Axis::Y) + end_shape(m, Axis::Y))) /
                                          6;
            }
    return start_shape.origin;
}

Vector3R gather_magnetic(const Field3d& field, const Vector3R& coord, const Domain& domain,
                         const BoundaryConditionHandler* bc) {
    const auto r = domain.to_cell_coordinates(coord);
    Vector3R result{0, 0, 0};
    for (int c = 0; c < 3; ++c) {
        const auto shift = magnetic_shift(c);
        const auto shifted = r - shift;
        Vector3I center;
        for (int d = 0; d < 3; ++d) center[d] = static_cast<int>(std::floor(shifted[d] + 0.5));
        const Stencil3D<3> shape = make_stencil_at<Shape2, 3>(shifted, center - Vector3I(1));
        int idx[3][3];
        for (int d = 0; d < 3; ++d)
            for (int n = 0; n < 3; ++n) idx[d][n] = node_owner(shape.origin[d] + n, d, domain, bc);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) {
                    const double w = shape(i, Axis::X) * shape(j, Axis::Y) * shape(k, Axis::Z);
                    result[c] += w * field(idx[0][i], idx[1][j], idx[2][k], c);
                }
    }
    return result;
}
}   // namespace

Vector3R gather_magnetic_point(const Field3d& magnetic, const Vector3R& position, const Domain& domain,
                               const BoundaryConditionHandler* bc) {
    return gather_magnetic(magnetic, position, domain, bc);
}

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
                    for (int c = 0; c < 3; ++c) current(ix, iy, iz, c) += value[c] * weights(n, m, k, c);
                }
            }
        }
    }
}

void deposit_charge(Field3d& density, const Vector3R& position, double charge_weight, const Domain& domain,
                    const BoundaryConditionHandler* bc) {
    const auto r = domain.to_cell_coordinates(wrap_position(position, domain, bc));
    Vector3I center;
    for (int d = 0; d < 3; ++d) center[d] = static_cast<int>(std::floor(r[d] + 0.5));
    const Stencil3D<3> shape = make_stencil_at<Shape2, 3>(r, center - Vector3I(1));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                const auto index = shape.origin + Vector3I(i, j, k);
                const double value = charge_weight * shape(i, Axis::X) * shape(j, Axis::Y) * shape(k, Axis::Z);
#pragma omp atomic update
                density(node_owner(index.x(), 0, domain, bc), node_owner(index.y(), 1, domain, bc),
                        node_owner(index.z(), 2, domain, bc), 0) += value;
            }
}

Particle advance_particle(const Particle& initial, const Field3d& electric, const Field3d& magnetic, Field3d& current,
                          double charge_mass, double charge_weight, double dt, const Domain& domain, double tolerance,
                          int max_iterations, const BoundaryConditionHandler* bc, const PathPolicy& policy,
                          OrbitBoundaryEvents* events, const Particle* guess) {
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
    struct Retry {
        double tau, position_error, velocity_residual, wall_fraction;
        int iterations;
        bool position_converged;
    };
    std::array<Retry, 60> retries;
    int retry_count = 0;
    OrbitFields last_fields;
    const auto vector_json = [](const Vector3R& v) { return nlohmann::json::array({v.x(), v.y(), v.z()}); };
    const auto particle_json = [&](const Particle& p) {
        return nlohmann::json{{"position", vector_json(p.coord)}, {"velocity", vector_json(p.velocity)}};
    };
    const auto fail = [&](const std::string& reason, double tau) -> void {
        nlohmann::json history = nlohmann::json::array();
        for (int i = 0; i < retry_count; ++i) {
            const auto& r = retries[i];
            history.push_back({{"tau", r.tau},
                               {"position_error", r.position_error},
                               {"velocity_residual", r.velocity_residual},
                               {"wall_fraction", r.wall_fraction},
                               {"iterations", r.iterations},
                               {"position_converged", r.position_converged}});
        }
        throw ParticleSolveFailure({{"reason", reason},
                                    {"initial", particle_json(initial)},
                                    {"substep_start", particle_json(start)},
                                    {"trial_end", particle_json(end)},
                                    {"dt", dt},
                                    {"remaining", remaining},
                                    {"tau", tau},
                                    {"substeps", substeps},
                                    {"tolerance", tolerance},
                                    {"max_iterations", max_iterations},
                                    {"charge_mass", charge_mass},
                                    {"charge_weight", charge_weight},
                                    {"last_gather_E", vector_json(last_fields.electric)},
                                    {"last_gather_B", vector_json(last_fields.magnetic)},
                                    {"retries", history}});
    };
    while (remaining > 0) {
        if (++substeps > 10000)
            fail("substep_limit", remaining);
        retry_count = 0;
        double tau = remaining;
        // Bracket the first wall hit between an interior endpoint (tau=0 at
        // first) and an exterior endpoint. Chord rescaling alone can stall
        // when the normal velocity approaches zero near a turning point.
        double wall_time_low = 0, wall_time_high = 0;
        bool converged = false;
        BoundaryCrossing hit;
        // A particle exactly on the wall and moving out interacts at tau=0.
        Vector3R clipped;
        hit = shape_particle_path(start.coord, start.coord + tau * start.velocity, domain, policy, clipped);
        if (!(hit.clipped && hit.t_hit == 0)) {
            hit = {};
            for (int retry = 0; retry < 60 && !converged; ++retry) {
                auto& trace = retries[retry];
                trace = {tau, -1, -1, -1, 0, false};
                retry_count = retry + 1;
                end = start;
                end.coord += tau * start.velocity;
                if (guess && substeps == 1 && retry == 0 && tau == dt) {
                    // The guess was wrapped through periodic faces: restore the orbit image.
                    Vector3R displacement = guess->coord - start.coord;
                    for (int d = 0; d < 3; ++d)
                        if (!bc || bc->is_periodic(d)) {
                            const double length = domain.num_cells(d) * domain.cell_size(d);
                            displacement[d] -= length * std::round(displacement[d] / length);
                        }
                    if (std::isfinite(displacement.length()) && std::isfinite(guess->velocity.length())) {
                        end.coord = start.coord + displacement;
                        end.velocity = guess->velocity;
                    }
                }
                OrbitAnderson acceleration;
                for (int it = 0; it < max_iterations; ++it) {
                    // Trial orbits can cross a wall. Gather only up to that wall;
                    // shorten tau below and solve the Lorentz equation again.
                    shape_particle_path(start.coord, end.coord, domain, policy, clipped);
                    const auto fields = gather_orbit(electric, magnetic, start.coord, clipped, domain, bc);
                    last_fields = fields;
                    trace.iterations = it + 1;
                    // CN image v_h = G(trial orbit).
                    const auto b = (0.5 * tau * charge_mass) * fields.magnetic;
                    const auto w = start.velocity + (0.5 * tau * charge_mass) * fields.electric;
                    const auto vh = (w + w.cross(b) + b * w.dot(b)) / (1 + b.dot(b));
                    const auto next_position = start.coord + tau * vh;
                    const auto next_velocity = 2.0 * vh - start.velocity;
                    const double error =
                        (next_position - end.coord).length() + tau * (next_velocity - end.velocity).length();
                    trace.position_error = error;
                    if (!std::isfinite(error))
                        fail("nonfinite_iterate", tau);
                    if (error <= tolerance * (1 + tau * vh.length())) {
                        end.coord = next_position;
                        end.velocity = next_velocity;
                        converged = true;
                        break;
                    }
                    const Vector3R trial = (end.coord - start.coord) / tau;
                    const Vector3R accelerated = acceleration.next(trial, vh);
                    end.coord = start.coord + tau * accelerated;
                    end.velocity = 2.0 * accelerated - start.velocity;
                }
                const auto crossing = shape_particle_path(start.coord, end.coord, domain, policy, clipped);
                trace.position_converged = converged;
                trace.wall_fraction = crossing.clipped ? crossing.t_hit : -1;
                if (converged && crossing.clipped) {
                    if ((end.coord - clipped).length() <= tolerance * (1 + end.coord.length())) {
                        end.coord = clipped;
                        hit = crossing;
                    } else {
                        const bool first_wall_trial = wall_time_high == 0;
                        wall_time_high = tau;
                        // The chord estimate is exact for free flight. After
                        // that first estimate, bisection contracts the bracket
                        // even for a strongly decelerating impact.
                        tau = first_wall_trial ? tau * crossing.t_hit
                                               : wall_time_low + 0.5 * (wall_time_high - wall_time_low);
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
                if (converged && !hit.clipped && wall_time_high > 0) {
                    wall_time_low = tau;
                    tau = wall_time_low + 0.5 * (wall_time_high - wall_time_low);
                    converged = false;
                    continue;
                }
                if (converged) {
                    const auto accepted = gather_orbit(electric, magnetic, start.coord, end.coord, domain, bc);
                    const auto vh = 0.5 * (start.velocity + end.velocity);
                    const auto residual = end.velocity - start.velocity -
                                          (tau * charge_mass) * (accepted.electric + vh.cross(accepted.magnetic));
                    last_fields = accepted;
                    trace.velocity_residual = residual.length();
                    converged = residual.length() <= tolerance * (1 + end.velocity.length());
                }
                if (converged)
                    break;
                hit = {};
                // A failed Lorentz solve changes the trial substep itself;
                // previous wall-time bounds no longer define this search.
                wall_time_low = wall_time_high = 0;
                tau *= 0.5;
            }
            if (!converged)
                fail("retry_limit", tau);
            if (!(tau > 0))
                fail("nonpositive_substep", tau);
            if (remaining - tau == remaining)
                fail("time_stagnation", tau);
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
            const int axis = (hit.face == Face::XMIN || hit.face == Face::XMAX)   ? X
                             : (hit.face == Face::YMIN || hit.face == Face::YMAX) ? Y
                                                                                  : Z;
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
