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
    double dt, qm, qw;
};

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

// Exact flow of H_Axis for a particle at cell + fraction, in cell units.
// The transverse position is fixed, so the path is cut only at the plane it crosses.
template <int Axis>
bool advance_axis(BucketFlow& flow, Particle& p, Vector3I& cell, Vector3R& fraction) {
    constexpr int A = (Axis + 1) % 3, B = (Axis + 2) % 3;
    const double h = flow.cell_size[Axis];
    const double coordinate = p.coord[Axis] + flow.dt * p.velocity[Axis];

    // End of the path relative to the starting cell: it must lie in this cell or in
    // a neighbor, as the bucket update requires. The negated test also rejects NaN.
    const double offset = coordinate / h - cell[Axis];
    if (!(offset >= -1 && offset < 2))
        return false;
    const int step = static_cast<int>(std::floor(offset));
    const double end_fraction = offset - step;

    // Path lengths, in cells, inside the starting cell and inside its neighbor.
    double length[2] = {end_fraction - fraction[Axis], 0};
    if (step > 0) {
        length[0] = 1 - fraction[Axis];
        length[1] = end_fraction;
    } else if (step < 0) {
        length[0] = -fraction[Axis];
        length[1] = end_fraction - 1;
    }

    const double wa[2] = {1 - fraction[A], fraction[A]};
    const double wb[2] = {1 - fraction[B], fraction[B]};
    Vector3I node = cell - flow.origin;
    for (int segment = 0; segment <= std::abs(step); ++segment) {
        const double path = length[segment] * h;

        // B_A is linear along A and constant across it; likewise B_B.
        Vector3I next_a = node, next_b = node;
        ++next_a[A];
        ++next_b[B];
        const double ba = wa[0] * flow.magnetic(node.x(), node.y(), node.z(), A) +
                          wa[1] * flow.magnetic(next_a.x(), next_a.y(), next_a.z(), A);
        const double bb = wb[0] * flow.magnetic(node.x(), node.y(), node.z(), B) +
                          wb[1] * flow.magnetic(next_b.x(), next_b.y(), next_b.z(), B);
        p.velocity[A] -= flow.qm * path * bb;
        p.velocity[B] += flow.qm * path * ba;

        for (int ia = 0; ia < 2; ++ia)
            for (int ib = 0; ib < 2; ++ib) {
                Vector3I edge = node;
                edge[A] += ia;
                edge[B] += ib;
                flow.current(edge.x(), edge.y(), edge.z(), Axis) += flow.qw * path * wa[ia] * wb[ib];
            }
        node[Axis] += step;
    }

    p.coord[Axis] = coordinate;
    cell[Axis] += step;
    fraction[Axis] = end_fraction;
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
                flow.dt = dt;
                flow.qm = qm;
                flow.qw = qw;
                load_patch(B, storage_origin, flow.magnetic);
                flow.current.zero();

                for (Particle& p : bucket) {
                    // The bucket fixes the cell: no search, and a patch access cannot leave the patch.
                    Vector3I cell = bucket_cell;
                    const Vector3R u = domain.to_cell_coordinates(p.coord);
                    Vector3R fraction(u.x() - cell.x(), u.y() - cell.y(), u.z() - cell.z());
                    // The subflows of one particle do not commute: B rotates the velocity in between.
                    if (reverse)
                        valid = valid && advance_axis<2>(flow, p, cell, fraction) &&
                                advance_axis<1>(flow, p, cell, fraction) && advance_axis<0>(flow, p, cell, fraction);
                    else
                        valid = valid && advance_axis<0>(flow, p, cell, fraction) &&
                                advance_axis<1>(flow, p, cell, fraction) && advance_axis<2>(flow, p, cell, fraction);
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
