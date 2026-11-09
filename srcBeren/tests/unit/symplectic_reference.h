#pragma once

#include "simulation/symplectic/whitney_pic.h"

namespace symplectic_reference {
inline int index(const Domain& domain, int i, int axis) {
    const int n = domain.num_cells(axis);
    return ((i % n) + n) % n + GHOST_CELLS;
}

inline double value(const Domain& domain, const Field3d& field, Vector3I node, int d) {
    return field(index(domain, node[0], 0), index(domain, node[1], 1), index(domain, node[2], 2), d);
}

inline Vector3R gather_reference(const Domain& domain, const Field3d& field, const Vector3R& position, bool magnetic) {
    const Vector3R u = domain.to_cell_coordinates(position);
    Vector3I cell;
    Vector3R f;
    for (int a = 0; a < 3; ++a) {
        cell[a] = static_cast<int>(std::floor(u[a]));
        f[a] = u[a] - cell[a];
    }
    Vector3R result(0, 0, 0);
    for (int d = 0; d < 3; ++d)
        for (int mask = 0; mask < 8; ++mask) {
            Vector3I node = cell;
            double w = 1;
            for (int a = 0; a < 3; ++a) {
                const int bit = (mask >> a) & 1;
                // E_d: constant in d, hat in other axes. B_d: hat in d,
                // constant in the two transverse axes. These are NOT CIC
                // weights centered independently at each Yee component.
                const bool linear = magnetic ? a == d : a != d;
                if (linear) {
                    node[a] += bit;
                    w *= bit ? f[a] : 1 - f[a];
                } else if (bit) {
                    w = 0;
                }
            }
            result[d] += w * value(domain, field, node, d);
        }
    return result;
}

inline void deposit_charge_reference(const Domain& domain, Field3d& rho, const Vector3R& position, double weight) {
    const Vector3R u = domain.to_cell_coordinates(position);
    Vector3I cell;
    Vector3R f;
    for (int a = 0; a < 3; ++a) {
        cell[a] = static_cast<int>(std::floor(u[a]));
        f[a] = u[a] - cell[a];
    }
    for (int mask = 0; mask < 8; ++mask) {
        Vector3I node = cell;
        double w = weight;
        for (int a = 0; a < 3; ++a) {
            const int bit = (mask >> a) & 1;
            node[a] = index(domain, cell[a] + bit, a);
            w *= bit ? f[a] : 1 - f[a];
        }
        rho(node[0], node[1], node[2], 0) += w;
    }
}

inline void kinetic_flow_reference(const Domain& domain, Particle& p, int axis, double dt, double qm, double qw,
                                   Field3d& E, const Field3d& B, Field3d& integral) {
    // v_axis is constant under H_axis. Split at EVERY crossed grid plane,
    // including periodic seams, so face-B and edge-current integrals are exact.
    const double h = domain.cell_size(axis);
    const double end = p.coord[axis] + dt * p.velocity[axis];
    double x = p.coord[axis];
    if (!std::isfinite(end))
        throw std::runtime_error("non-finite symplectic trajectory");
    while (x != end) {
        const bool forward = end > x;
        double plane = forward ? (std::floor(x / h) + 1) * h : (std::ceil(x / h) - 1) * h;
        // A previously reached plane may round to the preceding cell in x/h.
        if (forward ? plane <= x : plane >= x)
            plane += forward ? h : -h;
        const double next = forward ? std::min(end, plane) : std::max(end, plane);
        if (next == x)
            throw std::runtime_error("symplectic trajectory lost floating-point resolution");
        Vector3R midpoint = p.coord;
        midpoint[axis] = x + 0.5 * (next - x);
        const double length = next - x;
        Vector3R displacement(0, 0, 0);
        displacement[axis] = length;
        p.velocity += qm * displacement.cross(gather_reference(domain, B, midpoint, true));
        const Vector3R u = domain.to_cell_coordinates(midpoint);
        Vector3I cell;
        Vector3R f;
        for (int a = 0; a < 3; ++a) {
            cell[a] = static_cast<int>(std::floor(u[a]));
            f[a] = u[a] - cell[a];
        }
        const int a = (axis + 1) % 3, b = (axis + 2) % 3;
        for (int ia = 0; ia < 2; ++ia)
            for (int ib = 0; ib < 2; ++ib) {
                Vector3I node = cell;
                node[a] += ia;
                node[b] += ib;
                const double work_current = qw * length * (ia ? f[a] : 1 - f[a]) * (ib ? f[b] : 1 - f[b]);
                const int i = index(domain, node[0], 0), j = index(domain, node[1], 1), k = index(domain, node[2], 2);
                E(i, j, k, axis) -= work_current;
                integral(i, j, k, axis) += work_current;
            }
        x = next;
    }
    const double extent = domain.geom.box_max[axis] - domain.geom.box_min[axis];
    p.coord[axis] = end - std::floor(end / extent) * extent;
    if (p.coord[axis] >= extent)
        p.coord[axis] = 0;
}
}   // namespace symplectic_reference
