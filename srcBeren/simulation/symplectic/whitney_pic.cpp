#include "whitney_pic.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <source_location>

#include "coupling/esirkepov.h"
#include "coupling/gather.h"
#include "runtime/profiler.h"

namespace symplectic {
namespace {
// A particle starts in its bucket cell and every subflow crosses at most one plane,
// so it stays within one cell of the bucket: four nodes per axis, starting one below.
constexpr int PATCH_SIZE = 4;
using Patch = CurrentBuffer<PATCH_SIZE>;

// What the particles of one bucket share while B is fixed.
struct BucketFlow {
    Patch magnetic;
    Patch current;
    Vector3I origin;   // logical node of patch node (0, 0, 0)
    Vector3R cell_size;
    Vector3R inverse_cell_size;
    double dt, qm, qw;
};

// Distances between neighboring patch nodes in the flat patch storage.
constexpr int PATCH_STRIDE[3] = {Patch::idx(1, 0, 0, 0), Patch::idx(0, 1, 0, 0), Patch::idx(0, 0, 1, 0)};

// The part of a patch that lies in storage. Only a particle lying exactly on the
// last stored plane addresses the nodes beyond it, and then with zero weight.
struct PatchRange {
    Vector3I begin, end;

    PatchRange(const Field3d& field, const Vector3I& storage_origin) {
        for (int d = 0; d < 3; ++d) {
            begin[d] = std::max(0, -storage_origin[d]);
            end[d] = std::min(PATCH_SIZE, field.sizes()[d] - storage_origin[d]);
        }
    }
};

void load_patch(const Field3d& field, const Vector3I& storage_origin, Patch& patch) {
    const PatchRange range(field, storage_origin);
    patch.zero();
    for (int x = range.begin.x(); x < range.end.x(); ++x)
        for (int y = range.begin.y(); y < range.end.y(); ++y)
            for (int z = range.begin.z(); z < range.end.z(); ++z)
                for (int d = 0; d < 3; ++d)
                    patch(x, y, z, d) =
                        field(storage_origin.x() + x, storage_origin.y() + y, storage_origin.z() + z, d);
}

// Patches of neighboring buckets overlap. Slow particles leave most of a patch
// untouched, so the atomic update is paid only for the edges that received current.
void flush_patch(const Patch& patch, const Vector3I& storage_origin, Field3d& field) {
    const PatchRange range(field, storage_origin);
    for (int x = range.begin.x(); x < range.end.x(); ++x)
        for (int y = range.begin.y(); y < range.end.y(); ++y)
            for (int z = range.begin.z(); z < range.end.z(); ++z)
                for (int d = 0; d < 3; ++d) {
                    const double value = patch(x, y, z, d);
                    if (value == 0)
                        continue;
#pragma omp atomic update
                    field(storage_origin.x() + x, storage_origin.y() + y, storage_origin.z() + z, d) += value;
                }
}

// Exact flow of H_axis for a particle at cell + fraction, in cell units; node is the
// flat patch index of the cell. The transverse position is fixed, so the path is cut
// only at the plane it crosses.
template <int axis>
bool advance_axis(BucketFlow& flow, Particle& p, Vector3I& cell, int& node, Vector3R& fraction) {
    constexpr int a = (axis + 1) % 3, b = (axis + 2) % 3;
    const double h = flow.cell_size[axis];
    const double coordinate = p.coord[axis] + flow.dt * p.velocity[axis];

    // End of the path relative to the starting cell: it must lie in this cell or in
    // a neighbor, as the bucket update requires. The negated test also rejects NaN.
    const double offset = coordinate * flow.inverse_cell_size[axis] - cell[axis];
    if (!(offset >= -1 && offset < 2))
        return false;
    // offset + 1 is not negative, so truncation is the floor.
    const int step = static_cast<int>(offset + 1) - 1;
    const double end_fraction = offset - step;

    // Path lengths inside the starting cell and inside its neighbor.
    double path[2] = {(end_fraction - fraction[axis]) * h, 0};
    if (step > 0) {
        path[0] = (1 - fraction[axis]) * h;
        path[1] = end_fraction * h;
    } else if (step < 0) {
        path[0] = -fraction[axis] * h;
        path[1] = (end_fraction - 1) * h;
    }

    const double fa = fraction[a], fb = fraction[b];
    const double* magnetic = flow.magnetic.data_.data();
    double* current = flow.current.data_.data();
    int segment_node = node;
    for (int segment = 0; segment <= std::abs(step); ++segment) {
        // B_a is linear along a and constant across it; likewise B_b.
        const double ba = (1 - fa) * magnetic[segment_node + a] + fa * magnetic[segment_node + PATCH_STRIDE[a] + a];
        const double bb = (1 - fb) * magnetic[segment_node + b] + fb * magnetic[segment_node + PATCH_STRIDE[b] + b];
        p.velocity[a] -= flow.qm * path[segment] * bb;
        p.velocity[b] += flow.qm * path[segment] * ba;

        const int edge = segment_node + axis;
        const double charge = flow.qw * path[segment];
        current[edge] += charge * (1 - fa) * (1 - fb);
        current[edge + PATCH_STRIDE[b]] += charge * (1 - fa) * fb;
        current[edge + PATCH_STRIDE[a]] += charge * fa * (1 - fb);
        current[edge + PATCH_STRIDE[a] + PATCH_STRIDE[b]] += charge * fa * fb;
        segment_node += step * PATCH_STRIDE[axis];
    }

    p.coord[axis] = coordinate;
    cell[axis] += step;
    node += step * PATCH_STRIDE[axis];
    fraction[axis] = end_fraction;
    return true;
}
}   // namespace

void kick_velocities(ParticlesArray& particles, const Field3d& E, double impulse) {
    timer::commonTimer timer(std::source_location::current().function_name(), -1, timer::MeasureUnit::byte);
    const Domain& domain = particles.get_domain();
    int64_t count = 0;

#pragma omp parallel for schedule(dynamic, particles.omp_granularity()) reduction(+ : count)
    for (int bucket = 0; bucket < particles.size(); ++bucket) {
        count += std::ssize(particles(bucket));
        for (Particle& p : particles(bucket)) {
            const Vector3R coord = domain.to_cell_coordinates(p.coord);
            const Stencil3D<2> stencil = make_linear_stencil_at(coord, checked_floor_index(coord));
            p.velocity += impulse * gather_whitney_electric(E, stencil);
        }
    }
    timer.flat.m = sizeof(Particle) * count;
}

bool kinetic_flow(ParticlesArray& particles, const Field3d& B, double dt, bool reverse, Field3d& current) {
    timer::commonTimer timer(std::source_location::current().function_name(), -1, timer::MeasureUnit::byte);
    const Domain& domain = particles.get_domain();
    const Vector3I size = domain.size();
    const double qm = particles.charge() / particles.mass();
    const double qw = particles.charge() * particles.mpw();
    const Vector3R inverse_cell_size = Vector3R(1.0) / domain.cell_size();
    int64_t count = 0;
    bool valid = true;

#pragma omp parallel for collapse(3) schedule(dynamic, particles.omp_granularity()) reduction(+ : count) \
    reduction(&& : valid)
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int k = 0; k < size.z(); ++k) {
                std::vector<Particle>& bucket = particles(i, j, k);
                if (bucket.empty())
                    continue;
                count += std::ssize(bucket);

                const Vector3I bucket_cell = to_node_index(Vector3I(i, j, k));
                const Vector3I storage_origin(i - 1, j - 1, k - 1);
                BucketFlow flow;
                flow.origin = to_node_index(storage_origin);
                flow.cell_size = domain.cell_size();
                flow.inverse_cell_size = inverse_cell_size;
                flow.dt = dt;
                flow.qm = qm;
                flow.qw = qw;
                load_patch(B, storage_origin, flow.magnetic);
                flow.current.zero();

                for (Particle& p : bucket) {
                    // The bucket fixes the cell: no search, and a patch access cannot leave the patch.
                    Vector3I cell = bucket_cell;
                    int node = Patch::idx(1, 1, 1, 0);
                    // Same conversion as the path ends: a particle at rest must have a zero path.
                    const Vector3R u(p.coord.x() * flow.inverse_cell_size.x(), p.coord.y() * flow.inverse_cell_size.y(),
                                     p.coord.z() * flow.inverse_cell_size.z());
                    Vector3R fraction(u.x() - cell.x(), u.y() - cell.y(), u.z() - cell.z());
                    // The subflows of one particle do not commute: B rotates the velocity in between.
                    if (reverse)
                        valid = valid && advance_axis<2>(flow, p, cell, node, fraction) &&
                                advance_axis<1>(flow, p, cell, node, fraction) &&
                                advance_axis<0>(flow, p, cell, node, fraction);
                    else
                        valid = valid && advance_axis<0>(flow, p, cell, node, fraction) &&
                                advance_axis<1>(flow, p, cell, node, fraction) &&
                                advance_axis<2>(flow, p, cell, node, fraction);
                }
                flush_patch(flow.current, storage_origin, current);
            }
    timer.flat.m = sizeof(Particle) * count;
    return valid;
}

double max_magnetic_divergence(const Domain& domain, const Field3d& B) {
    RECORD_TIMER;
    const Vector3I size = domain.size();
    const Vector3R h = domain.cell_size();
    double result = 0;

    // B_d(i) lies on the lower d-face of cell i; the last stored plane has no upper neighbor.
#pragma omp parallel for collapse(3) reduction(max : result)
    for (int i = 0; i < size.x() - 1; ++i)
        for (int j = 0; j < size.y() - 1; ++j)
            for (int k = 0; k < size.z() - 1; ++k) {
                const double divergence = (B(i + 1, j, k, 0) - B(i, j, k, 0)) / h.x() +
                                          (B(i, j + 1, k, 1) - B(i, j, k, 1)) / h.y() +
                                          (B(i, j, k + 1, 2) - B(i, j, k, 2)) / h.z();
                result = std::max(result, std::abs(divergence));
            }
    return result;
}
}   // namespace symplectic
