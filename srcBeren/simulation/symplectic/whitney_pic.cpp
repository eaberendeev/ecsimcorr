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

// Particles advanced together; 4 was faster than 8 and 16 on a Haswell desktop.
constexpr int BATCH_WIDTH = 4;

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

// Particles of one bucket advanced together, one lane per particle. Every lane runs
// the same code with the same index arithmetic. A padding lane stands still, so its
// paths, and with them all its contributions, are exactly zero.
template <int width>
struct ParticleBatch {
    double coord[3][width];
    double velocity[3][width];
    double fraction[3][width];   // position inside the cell, in cell units
    int node[3][width];          // patch node of the cell
};

// Distances between neighboring patch nodes in the flat patch storage.
constexpr int PATCH_STRIDE[3] = {Patch::idx(1, 0, 0, 0), Patch::idx(0, 1, 0, 0), Patch::idx(0, 0, 1, 0)};

// Rotation of the transverse velocity and current of one path segment starting at patch node `node`.
template <int axis, int width>
void add_segment(BucketFlow& flow, ParticleBatch<width>& batch, int lane, int node, double rotation, double charge) {
    constexpr int a = (axis + 1) % 3, b = (axis + 2) % 3;
    const double fa = batch.fraction[a][lane], fb = batch.fraction[b][lane];
    const double* magnetic = flow.magnetic.data_.data();
    double* current = flow.current.data_.data();

    // B_a is linear along a and constant across it; likewise B_b.
    const double ba = (1 - fa) * magnetic[node + a] + fa * magnetic[node + PATCH_STRIDE[a] + a];
    const double bb = (1 - fb) * magnetic[node + b] + fb * magnetic[node + PATCH_STRIDE[b] + b];
    batch.velocity[a][lane] -= rotation * bb;
    batch.velocity[b][lane] += rotation * ba;

    const int edge = node + axis;
    current[edge] += charge * (1 - fa) * (1 - fb);
    current[edge + PATCH_STRIDE[b]] += charge * (1 - fa) * fb;
    current[edge + PATCH_STRIDE[a]] += charge * fa * (1 - fb);
    current[edge + PATCH_STRIDE[a] + PATCH_STRIDE[b]] += charge * fa * fb;
}

// Exact flow of H_axis for every lane. The transverse position is fixed, so a path
// is cut only at the plane it crosses: one segment in the starting cell and, for the
// lanes that cross, one in the neighbor.
template <int axis, int width>
bool advance_axis(BucketFlow& flow, ParticleBatch<width>& batch) {
    // Locals: the batch stores must not look as if they could change the flow.
    const double h = flow.cell_size[axis], inverse_h = 1 / h, dt = flow.dt;
    const double origin = flow.origin[axis];
    double coordinate[width], offset[width];
    bool valid = true;

    // End of the path relative to the starting cell: it must lie in this cell or in
    // a neighbor, as the bucket update requires. The comparisons also reject NaN.
    for (int lane = 0; lane < width; ++lane) {
        coordinate[lane] = batch.coord[axis][lane] + dt * batch.velocity[axis][lane];
        offset[lane] = coordinate[lane] * inverse_h - (origin + batch.node[axis][lane]);
        valid &= (offset[lane] >= -1) & (offset[lane] < 2);
    }
    if (!valid)
        return false;

    int step[width], start[width];
    double path[2][width];
    int crossing[width];
    int crossings = 0;
    for (int lane = 0; lane < width; ++lane) {
        // offset + 1 is not negative, so truncation is the floor.
        step[lane] = static_cast<int>(offset[lane] + 1) - 1;
        const double end = offset[lane] - step[lane];
        const double inside = step[lane] == 0 ? end : (step[lane] > 0 ? 1.0 : 0.0);
        const double beyond = step[lane] == 0 ? 0.0 : (step[lane] > 0 ? end : end - 1);
        path[0][lane] = (inside - batch.fraction[axis][lane]) * h;
        path[1][lane] = beyond * h;
        batch.coord[axis][lane] = coordinate[lane];
        batch.fraction[axis][lane] = end;

        start[lane] = PATCH_STRIDE[0] * batch.node[0][lane] + PATCH_STRIDE[1] * batch.node[1][lane] +
                      PATCH_STRIDE[2] * batch.node[2][lane];
        batch.node[axis][lane] += step[lane];
        crossing[crossings] = lane;
        crossings += step[lane] != 0;
    }
    const double qm = flow.qm, qw = flow.qw;

    // The segment in the starting cell, for every lane.
    for (int lane = 0; lane < width; ++lane)
        add_segment<axis>(flow, batch, lane, start[lane], qm * path[0][lane], qw * path[0][lane]);
    // The segment in the neighbor, only for the lanes that crossed a plane.
    for (int k = 0; k < crossings; ++k) {
        const int lane = crossing[k];
        add_segment<axis>(flow, batch, lane, start[lane] + step[lane] * PATCH_STRIDE[axis], qm * path[1][lane],
                          qw * path[1][lane]);
    }
    return true;
}

// Lanes past the bucket end repeat its first particle at rest.
template <int width>
void load_batch(const Domain& domain, const std::vector<Particle>& bucket, int first, const Vector3I& bucket_cell,
                ParticleBatch<width>& batch) {
    const int count = std::min<int>(width, std::ssize(bucket) - first);
    for (int lane = 0; lane < width; ++lane) {
        const Particle& p = bucket[first + (lane < count ? lane : 0)];
        for (int d = 0; d < 3; ++d) {
            const double u_d = p.coord[d] * (1 / domain.cell_size(d));
            batch.coord[d][lane] = p.coord[d];
            batch.velocity[d][lane] = lane < count ? p.velocity[d] : 0.0;
            batch.fraction[d][lane] = u_d - bucket_cell[d];
            batch.node[d][lane] = 1;
        }
    }
}

template <int width>
void store_batch(const ParticleBatch<width>& batch, int first, std::vector<Particle>& bucket) {
    const int count = std::min<int>(width, std::ssize(bucket) - first);
    for (int lane = 0; lane < count; ++lane)
        for (int d = 0; d < 3; ++d) {
            bucket[first + lane].coord[d] = batch.coord[d][lane];
            bucket[first + lane].velocity[d] = batch.velocity[d][lane];
        }
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

                // The bucket fixes the cell: no search, and a patch access cannot leave the patch.
                for (int first = 0; first < std::ssize(bucket); first += BATCH_WIDTH) {
                    ParticleBatch<BATCH_WIDTH> batch;
                    load_batch(domain, bucket, first, bucket_cell, batch);
                    // The subflows of one particle do not commute: B rotates the velocity in between.
                    if (reverse)
                        valid = valid && advance_axis<2>(flow, batch) && advance_axis<1>(flow, batch) &&
                                advance_axis<0>(flow, batch);
                    else
                        valid = valid && advance_axis<0>(flow, batch) && advance_axis<1>(flow, batch) &&
                                advance_axis<2>(flow, batch);
                    store_batch(batch, first, bucket);
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
