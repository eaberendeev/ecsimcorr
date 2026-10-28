#include "whitney_pic.h"

#include <cmath>
#include <stdexcept>

namespace symplectic {
int WhitneyPIC::index(int i, int axis) const {
    const int n = domain_.num_cells(axis);
    return ((i % n) + n) % n + domain_.grid.ghost_cells();
}
double WhitneyPIC::value(const Field3d& field, Vector3I node, int component) const {
    return field(index(node[0], 0), index(node[1], 1), index(node[2], 2), component);
}
void WhitneyPIC::sync(Field3d& field) const {
    const int g = domain_.grid.ghost_cells();
    const auto size = domain_.size();
    for (int i = 0; i < size[0]; ++i)
        for (int j = 0; j < size[1]; ++j)
            for (int k = 0; k < size[2]; ++k)
                for (int d = 0; d < field.nd(); ++d) field(i, j, k, d) = value(field, Vector3I(i - g, j - g, k - g), d);
}
Vector3R WhitneyPIC::gather(const Field3d& field, const Vector3R& position, bool magnetic) const {
    const auto u = domain_.to_cell_coordinates(position);
    Vector3I cell;
    Vector3R f;
    for (int a = 0; a < 3; ++a) {
        cell[a] = static_cast<int>(std::floor(u[a]));
        f[a] = u[a] - cell[a];
    }
    Vector3R result(0, 0, 0);
    for (int d = 0; d < 3; ++d)
        for (int mask = 0; mask < 8; ++mask) {
            auto node = cell;
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
            result[d] += w * value(field, node, d);
        }
    return result;
}
void WhitneyPIC::deposit_charge(Field3d& rho, const Vector3R& position, double weight) const {
    const auto u = domain_.to_cell_coordinates(position);
    Vector3I cell;
    Vector3R f;
    for (int a = 0; a < 3; ++a) {
        cell[a] = static_cast<int>(std::floor(u[a]));
        f[a] = u[a] - cell[a];
    }
    for (int mask = 0; mask < 8; ++mask) {
        auto node = cell;
        double w = weight;
        for (int a = 0; a < 3; ++a) {
            const int bit = (mask >> a) & 1;
            node[a] = index(cell[a] + bit, a);
            w *= bit ? f[a] : 1 - f[a];
        }
        rho(node[0], node[1], node[2], 0) += w;
    }
}
void WhitneyPIC::electric_flow(Field3d& B, const Field3d& E, double dt) const {
    const int g = domain_.grid.ghost_cells();
    const auto n = domain_.num_cells();
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const int a = (d + 1) % 3, b = (d + 2) % 3;
                    Vector3I node(i, j, k), na = node, nb = node;
                    ++na[a];
                    ++nb[b];
                    const double curl = (value(E, na, b) - value(E, node, b)) / domain_.cell_size(a) -
                                        (value(E, nb, a) - value(E, node, a)) / domain_.cell_size(b);
                    B(i + g, j + g, k + g, d) -= dt * curl;
                }
}
void WhitneyPIC::magnetic_flow(Field3d& E, const Field3d& B, double dt) const {
    const int g = domain_.grid.ghost_cells();
    const auto n = domain_.num_cells();
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const int a = (d + 1) % 3, b = (d + 2) % 3;
                    Vector3I node(i, j, k), na = node, nb = node;
                    --na[a];
                    --nb[b];
                    const double curl = (value(B, node, b) - value(B, na, b)) / domain_.cell_size(a) -
                                        (value(B, node, a) - value(B, nb, a)) / domain_.cell_size(b);
                    E(i + g, j + g, k + g, d) += dt * curl;
                }
}
void WhitneyPIC::kinetic_flow(Particle& p, int axis, double dt, double qm, double qw, Field3d& E, const Field3d& B,
                              Field3d& integral) const {
    // v_axis is constant under H_axis. Split at EVERY crossed grid plane,
    // including periodic seams, so face-B and edge-current integrals are exact.
    const double h = domain_.cell_size(axis);
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
        auto midpoint = p.coord;
        midpoint[axis] = x + 0.5 * (next - x);
        const double length = next - x;
        Vector3R displacement(0, 0, 0);
        displacement[axis] = length;
        p.velocity += qm * displacement.cross(gather(B, midpoint, true));
        const auto u = domain_.to_cell_coordinates(midpoint);
        Vector3I cell;
        Vector3R f;
        for (int a = 0; a < 3; ++a) {
            cell[a] = static_cast<int>(std::floor(u[a]));
            f[a] = u[a] - cell[a];
        }
        const int a = (axis + 1) % 3, b = (axis + 2) % 3;
        for (int ia = 0; ia < 2; ++ia)
            for (int ib = 0; ib < 2; ++ib) {
                auto node = cell;
                node[a] += ia;
                node[b] += ib;
                const double work_current = qw * length * (ia ? f[a] : 1 - f[a]) * (ib ? f[b] : 1 - f[b]);
                const int i = index(node[0], 0), j = index(node[1], 1), k = index(node[2], 2);
                E(i, j, k, axis) -= work_current;
                integral(i, j, k, axis) += work_current;
            }
        x = next;
    }
    const double extent = domain_.geom.box_max[axis] - domain_.geom.box_min[axis];
    p.coord[axis] = end - std::floor(end / extent) * extent;
    if (p.coord[axis] >= extent)
        p.coord[axis] = 0;
}
double WhitneyPIC::divergence(const Field3d& field, int i, int j, int k, bool magnetic) const {
    const Vector3I node(i, j, k);
    double result = 0;
    for (int d = 0; d < 3; ++d) {
        auto other = node;
        other[d] += magnetic ? 1 : -1;
        result += (magnetic ? value(field, other, d) - value(field, node, d)
                            : value(field, node, d) - value(field, other, d)) /
                  domain_.cell_size(d);
    }
    return result;
}
double WhitneyPIC::field_energy(const Field3d& field) const {
    double sum = 0;
    const auto n = domain_.num_cells();
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const double v = value(field, Vector3I(i, j, k), d);
                    sum += 0.5 * v * v;
                }
    return sum;
}
}   // namespace symplectic
