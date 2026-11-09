#include "whitney_pic.h"

#include <cmath>
#include <stdexcept>

#include "runtime/profiler.h"

namespace symplectic {
namespace {
struct ParticleStencil {
    Vector3I cell, low, high;
    Vector3R fraction;

    ParticleStencil(const Domain& domain, const Vector3R& position) {
        const Vector3R u = domain.to_cell_coordinates(position);
        cell = checked_floor_index(u);
        fraction = u - Vector3R(cell[0], cell[1], cell[2]);
        for (int d = 0; d < 3; ++d) {
            low[d] = domain.grid.periodic_storage_index(cell[d], d);
            high[d] = domain.grid.periodic_storage_index(cell[d] + 1, d);
        }
    }

    void update_axis(const Domain& domain, const Vector3R& position, int axis) {
        const double u = position[axis] / domain.cell_size(axis);
        cell[axis] = checked_floor_index(u);
        fraction[axis] = u - cell[axis];
        low[axis] = domain.grid.periodic_storage_index(cell[axis], axis);
        high[axis] = domain.grid.periodic_storage_index(cell[axis] + 1, axis);
    }
};

// One bucket shares a 4x4x4 patch. The second half-step may leave this patch
// because particles are rebinned only at the end of the complete step.
// Current stays local until every particle in the bucket has been advanced.
struct CellPatch {
    Vector3I start, cells, size;
    double magnetic[64][3];
    double current[64][3] = {};

    CellPatch(const Domain& domain, int bucket, const Field3d& field) : cells(domain.num_cells()) {
        for (int d = 0; d < 3; ++d) {
            start[d] = domain.grid.periodic_storage_index(domain.pos_sind(bucket, d) - GHOST_CELLS - 1, d);
            size[d] = std::min(4, cells[d]);
        }
        for (int i = 0; i < size[0]; ++i)
            for (int j = 0; j < size[1]; ++j)
                for (int k = 0; k < size[2]; ++k) {
                    const Vector3I node = storage_node(i, j, k);
                    for (int d = 0; d < 3; ++d) magnetic[(i * 4 + j) * 4 + k][d] = field(node[0], node[1], node[2], d);
                }
    }

    Vector3I storage_node(int i, int j, int k) const {
        return {GHOST_CELLS + (start[0] - GHOST_CELLS + i) % cells[0],
                GHOST_CELLS + (start[1] - GHOST_CELLS + j) % cells[1],
                GHOST_CELLS + (start[2] - GHOST_CELLS + k) % cells[2]};
    }

    int local_index(int node, int axis) const {
        const int local = node - start[axis];
        return local < 0 ? local + cells[axis] : local;
    }

    void flush(Field3d& integral) const {
        for (int i = 0; i < size[0]; ++i)
            for (int j = 0; j < size[1]; ++j)
                for (int k = 0; k < size[2]; ++k) {
                    const Vector3I node = storage_node(i, j, k);
                    for (int d = 0; d < 3; ++d)
                        integral(node[0], node[1], node[2], d) += current[(i * 4 + j) * 4 + k][d];
                }
    }
};

// Each axis subflow crosses at most one grid plane.
double axis_displacement(const Domain& domain, const Particle& particle, int axis, double dt) {
    const double displacement = dt * particle.velocity[axis];
    if (!std::isfinite(displacement) || !std::isfinite(particle.coord[axis] + displacement))
        throw std::runtime_error("non-finite symplectic trajectory");
    if (std::abs(displacement) > domain.cell_size(axis))
        throw std::runtime_error("symplectic particle displacement exceeds one cell; reduce Dt");
    return displacement;
}

double wrap_coordinate(double coordinate, double extent) {
    if (coordinate < 0)
        coordinate += extent;
    else if (coordinate >= extent)
        coordinate -= extent;
    return coordinate >= extent ? 0 : coordinate;
}

}   // namespace

static void advance_axis_impl(const Domain& domain, Particle& p, int axis, double dt, double qm, double qw,
                              const Field3d& B, Field3d& integral, CellPatch* patch, ParticleStencil& stencil);

static double value(const Domain& domain, const Field3d& field, Vector3I node, int component) {
    return field(domain.grid.periodic_storage_index(node[0], 0), domain.grid.periodic_storage_index(node[1], 1),
                 domain.grid.periodic_storage_index(node[2], 2), component);
}

void fill_periodic_ghosts(const Domain& domain, Field3d& field) {
    RECORD_TIMER;
    const int g = GHOST_CELLS;
    const Vector3I size = domain.size();
    const Vector3I n = domain.num_cells();
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < size[0]; ++i)
        for (int j = 0; j < size[1]; ++j)
            for (int k = 0; k < size[2]; ++k) {
                if (i >= g && i < g + n[0] && j >= g && j < g + n[1] && k >= g && k < g + n[2])
                    continue;
                for (int d = 0; d < field.nd(); ++d)
                    field(i, j, k, d) = value(domain, field, Vector3I(i - g, j - g, k - g), d);
            }
}

template <bool isMagnetic>
Vector3R gather(const Domain& domain, const Field3d& field, const Vector3R& position) {
    const ParticleStencil stencil(domain, position);
    const Vector3I& low = stencil.low;
    const Vector3I& high = stencil.high;
    const Vector3R& fraction = stencil.fraction;
    Vector3R result(0, 0, 0);
    for (int d = 0; d < 3; ++d) {
        if constexpr (isMagnetic) {
            Vector3I upper = low;
            upper[d] = high[d];
            result[d] = (1 - fraction[d]) * field(low[0], low[1], low[2], d) +
                        fraction[d] * field(upper[0], upper[1], upper[2], d);
        } else {
            const int a = (d + 1) % 3, b = (d + 2) % 3;
            for (int ia = 0; ia < 2; ++ia)
                for (int ib = 0; ib < 2; ++ib) {
                    Vector3I node = low;
                    if (ia)
                        node[a] = high[a];
                    if (ib)
                        node[b] = high[b];
                    const double w = (ia ? fraction[a] : 1 - fraction[a]) * (ib ? fraction[b] : 1 - fraction[b]);
                    result[d] += w * field(node[0], node[1], node[2], d);
                }
        }
    }
    return result;
}

void deposit_charge(const Domain& domain, Field3d& rho, const Vector3R& position, double weight) {
    const ParticleStencil stencil(domain, position);
    const Vector3I& low = stencil.low;
    const Vector3I& high = stencil.high;
    const Vector3R& fraction = stencil.fraction;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                const double wx = i ? fraction[0] : 1 - fraction[0];
                const double wy = j ? fraction[1] : 1 - fraction[1];
                const double wz = k ? fraction[2] : 1 - fraction[2];
                rho(i ? high[0] : low[0], j ? high[1] : low[1], k ? high[2] : low[2], 0) += weight * wx * wy * wz;
            }
}

void advance_axis(const Domain& domain, Particle& p, int axis, double dt, double qm, double qw, const Field3d& B,
                  Field3d& integral) {
    if (qm == 0 && qw == 0) {
        const double displacement = axis_displacement(domain, p, axis, dt);
        p.coord[axis] = wrap_coordinate(p.coord[axis] + displacement, domain.cell_size(axis) * domain.num_cells(axis));
        return;
    }
    ParticleStencil stencil(domain, p.coord);
    advance_axis_impl(domain, p, axis, dt, qm, qw, B, integral, nullptr, stencil);
}

static void advance_axis_impl(const Domain& domain, Particle& p, int axis, double dt, double qm, double qw,
                              const Field3d& B, Field3d& integral, CellPatch* patch, ParticleStencil& stencil) {
    const double h = domain.cell_size(axis);
    const double displacement = axis_displacement(domain, p, axis, dt);
    const double end = p.coord[axis] + displacement;
    const double extent = domain.cell_size(axis) * domain.num_cells(axis);
    Vector3I& low = stencil.low;
    const Vector3I& high = stencil.high;
    const Vector3R& fraction = stencil.fraction;
    const int a = (axis + 1) % 3, b = (axis + 2) % 3;
    const int stride[3] = {16, 4, 1};
    const int local_a[2] = {patch ? patch->local_index(low[a], a) : 0, patch ? patch->local_index(high[a], a) : 0};
    const int local_b[2] = {patch ? patch->local_index(low[b], b) : 0, patch ? patch->local_index(high[b], b) : 0};
    const bool transverse_cached = patch && local_a[0] < 4 && local_a[1] < 4 && local_b[0] < 4 && local_b[1] < 4;
    const int offset_a = (local_a[1] - local_a[0]) * stride[a];
    const int offset_b = (local_b[1] - local_b[0]) * stride[b];
    const double start = p.coord[axis];
    if (start == end)
        return;

    const bool forward = displacement > 0;
    int along = stencil.cell[axis];
    double plane;
    if (forward) {
        plane = (along + 1) * h;
        if (plane <= start) {
            ++along;
            plane = (along + 1) * h;
        }
    } else {
        plane = along * h;
        if (plane >= start) {
            --along;
            plane = along * h;
        }
    }
    const bool crosses = forward ? end > plane : end < plane;
    const double first_end = crosses ? plane : end;
    const double lengths[2] = {first_end - start, end - first_end};
    // The transverse stencil stays fixed. Integrate over the starting cell and,
    // only when crossed, its neighbor; no cell search is needed within a segment.
    for (int segment = 0; segment < (crosses ? 2 : 1); ++segment) {
        const double length = lengths[segment];
        low[axis] = domain.grid.periodic_storage_index(along, axis);
        Vector3I upper_a = low, upper_b = low;
        upper_a[a] = high[a];
        upper_b[b] = high[b];
        const int local_axis = patch ? patch->local_index(low[axis], axis) : 0;
        const bool cached = transverse_cached && local_axis < 4;
        const int base = local_axis * stride[axis] + local_a[0] * stride[a] + local_b[0] * stride[b];
        double ba, bb;
        if (cached) {
            ba = (1 - fraction[a]) * patch->magnetic[base][a] + fraction[a] * patch->magnetic[base + offset_a][a];
            bb = (1 - fraction[b]) * patch->magnetic[base][b] + fraction[b] * patch->magnetic[base + offset_b][b];
        } else {
            ba = (1 - fraction[a]) * B(low[0], low[1], low[2], a) +
                 fraction[a] * B(upper_a[0], upper_a[1], upper_a[2], a);
            bb = (1 - fraction[b]) * B(low[0], low[1], low[2], b) +
                 fraction[b] * B(upper_b[0], upper_b[1], upper_b[2], b);
        }
        p.velocity[a] -= qm * (length * bb);
        p.velocity[b] += qm * (length * ba);
        for (int ia = 0; ia < 2; ++ia)
            for (int ib = 0; ib < 2; ++ib) {
                Vector3I node = low;
                if (ia)
                    node[a] = high[a];
                if (ib)
                    node[b] = high[b];
                const double current =
                    qw * length * (ia ? fraction[a] : 1 - fraction[a]) * (ib ? fraction[b] : 1 - fraction[b]);
                if (cached)
                    patch->current[base + ia * offset_a + ib * offset_b][axis] += current;
                else
                    integral(node[0], node[1], node[2], axis) += current;
            }
        along += forward ? 1 : -1;
    }
    p.coord[axis] = wrap_coordinate(end, extent);
    stencil.update_axis(domain, p.coord, axis);
}

void drift_neutral(const Domain& domain, Particle& particle, double dt) {
    Vector3R position = particle.coord;
    for (int axis = 0; axis < 3; ++axis) {
        const double displacement = axis_displacement(domain, particle, axis, dt);
        position[axis] =
            wrap_coordinate(position[axis] + displacement, domain.cell_size(axis) * domain.num_cells(axis));
    }
    particle.coord = position;
}

void advance_bucket(const Domain& domain, std::vector<Particle>& particles, int bucket, bool reverse, double dt,
                    double qm, double qw, const Field3d& magnetic, Field3d& integral) {
    if (particles.empty())
        return;

    if (qm == 0 && qw == 0) {
        for (Particle& particle : particles) drift_neutral(domain, particle, dt);
        return;
    }

    CellPatch patch(domain, bucket, magnetic);
    for (Particle& particle : particles) {
        ParticleStencil stencil(domain, particle.coord);
        // Different particles commute while B is fixed; axes of one particle do not.
        for (int step = 0; step < 3; ++step) {
            const int axis = reverse ? 2 - step : step;
            advance_axis_impl(domain, particle, axis, dt, qm, qw, magnetic, integral, &patch, stencil);
        }
    }
    patch.flush(integral);
}

void electric_flow(const Domain& domain, Field3d& B, const Field3d& E, double dt) {
    RECORD_TIMER;
    const int g = GHOST_CELLS;
    const Vector3I n = domain.num_cells();
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const int a = (d + 1) % 3, b = (d + 2) % 3;
                    Vector3I node(i, j, k), na = node, nb = node;
                    ++na[a];
                    ++nb[b];
                    const double curl = (value(domain, E, na, b) - value(domain, E, node, b)) / domain.cell_size(a) -
                                        (value(domain, E, nb, a) - value(domain, E, node, a)) / domain.cell_size(b);
                    B(i + g, j + g, k + g, d) -= dt * curl;
                }
}

void magnetic_flow(const Domain& domain, Field3d& E, const Field3d& B, double dt) {
    RECORD_TIMER;
    const int g = GHOST_CELLS;
    const Vector3I n = domain.num_cells();
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const int a = (d + 1) % 3, b = (d + 2) % 3;
                    Vector3I node(i, j, k), na = node, nb = node;
                    --na[a];
                    --nb[b];
                    const double curl = (value(domain, B, node, b) - value(domain, B, na, b)) / domain.cell_size(a) -
                                        (value(domain, B, node, a) - value(domain, B, nb, a)) / domain.cell_size(b);
                    E(i + g, j + g, k + g, d) += dt * curl;
                }
}

double divergence(const Domain& domain, const Field3d& field, int i, int j, int k, bool isMagnetic) {
    const Vector3I node(i, j, k);
    double result = 0;
    if (isMagnetic) {
        for (int d = 0; d < 3; ++d) {
            Vector3I next = node;
            ++next[d];
            result += (value(domain, field, next, d) - value(domain, field, node, d)) / domain.cell_size(d);
        }
    } else {
        for (int d = 0; d < 3; ++d) {
            Vector3I previous = node;
            --previous[d];
            result += (value(domain, field, node, d) - value(domain, field, previous, d)) / domain.cell_size(d);
        }
    }
    return result;
}

double field_energy(const Domain& domain, const Field3d& field) {
    RECORD_TIMER;
    double sum = 0;
    const Vector3I n = domain.num_cells();
    std::vector<double> slabs(n[0], 0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const double v = value(domain, field, Vector3I(i, j, k), d);
                    slabs[i] += 0.5 * v * v;
                }
    for (double slab : slabs) sum += slab;
    return sum;
}

template Vector3R gather<false>(const Domain&, const Field3d&, const Vector3R&);
template Vector3R gather<true>(const Domain&, const Field3d&, const Vector3R&);
}   // namespace symplectic
