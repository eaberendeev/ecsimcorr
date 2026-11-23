#pragma once

#include <omp.h>

#include <cmath>
#include <vector>

#include "coupling/gather.h"
#include "particles/particle_array.h"
#include "runtime/profiler.h"

// Missing magnetization of the straight step path. The path x^n -> x^{n+1} =
// x^n + dt v_mid is the Crank-Nicolson chord of the Boris rotation: a gyrating
// particle deposits on it only the magnetic moment mu/(1 + theta^2), theta =
// |q/m| |B| dt/2 (the two half-paths of ECSIM-Corr enclose the full one).
// Each particle adds at x^{n+1/2}, with the weights of the linear Yee
// magnetic gather,
//   dM_p = -mpw (mu - mu_eff) b,  mu - mu_eff = m u^2/(2|B|) theta^2/(1 + theta^2),
// u the gyration speed (v_perp - v_E at x^{n+1/2}, v^n); the field equations get
// J_M = curl_B dM, and the particle's gyration energy changes by
// -dM_p . (B^{n+1} - B^n) with the same weights. Since curl_B = curl_E^T and
// B^{n+1} - B^n = -dt curl_E E_mean, the field loses what the particles gain
// (the guiding-centre theorem dK = dt E.J - M.dB): the energy balance is kept.
// At Omega dt << 1, dM -> 0; at B = 0 it vanishes.
namespace pathMagnetization {

// Moment magnitudes mpw (mu - mu_eff) per particle, in the order of the cells.
using Moments = Array3D<std::vector<double>>;

struct Exchange {
    double target = 0, given = 0;
};

// Weights of the linear Yee magnetic gather at a position, visit(i, j, k,
// component, weight) on logical nodes, halo included: a deposit is folded like
// a current.
template <typename Visit>
inline void visit_magnetic_weights(const Vector3R& position, const Domain& domain, Visit&& visit) {
    const Vector3R coord = domain.to_cell_coordinates(position);
    const Vector3I origin = checked_floor_index(coord);
    const Stencil3D<2> no = make_linear_stencil_at(coord, origin);
    const Stencil3D<2> sh =
        make_linear_stencil_at(coord, linear_half_shifted_origin(coord, origin), Vector3R::Constant(0.5));
    for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x) {
                const int ix = no.node(x, Axis::X), iy = no.node(y, Axis::Y), iz = no.node(z, Axis::Z);
                const int hx = sh.node(x, Axis::X), hy = sh.node(y, Axis::Y), hz = sh.node(z, Axis::Z);
                visit(ix, hy, hz, 0, no(x, Axis::X) * sh(y, Axis::Y) * sh(z, Axis::Z));
                visit(hx, iy, hz, 1, sh(x, Axis::X) * no(y, Axis::Y) * sh(z, Axis::Z));
                visit(hx, hy, iz, 2, sh(x, Axis::X) * sh(y, Axis::Y) * no(z, Axis::Z));
            }
}

// Unit field direction at a position (zero where B = 0).
inline Vector3R field_direction(const Field3d& magnetic, const Vector3R& position, const Domain& domain) {
    const Vector3R coord = domain.to_cell_coordinates(position);
    const Vector3I origin = checked_floor_index(coord);
    const Stencil3D<2> no = make_linear_stencil_at(coord, origin);
    const Stencil3D<2> sh =
        make_linear_stencil_at(coord, linear_half_shifted_origin(coord, origin), Vector3R::Constant(0.5));
    const Vector3R field = interpolateB(magnetic, no, sh);
    const double strength = field.norm();
    return strength > 0 ? Vector3R(field / strength) : Vector3R(Vector3R::Zero());
}

// Particles at x^{n+1/2} with v^n: their missing moments into `moments` and
// dM into `moment` (B layout).
inline void deposit(const ParticlesArray& sp, const Field3d& magnetic, const Field3d& electric, double dt,
                    Moments& moments, Field3d& moment) {
    RECORD_TIMER;
    const Domain& domain = sp.get_domain();
    const double qm = sp.charge() / sp.mass(), mw = sp.mass() * sp.mpw();
    moments.resize(sp.cells().size());
    const int threads = omp_get_max_threads();
    std::vector<Field3d> partial(threads, Field3d(moment.sizes(), 3));
#pragma omp parallel num_threads(threads)
    {
        Field3d& grid = partial[omp_get_thread_num()];
        grid.setZero();
        const GridNodeView<Field3d> view(grid);
#pragma omp for schedule(dynamic, sp.omp_granularity())
        for (int cell = 0; cell < sp.size(); ++cell) {
            const std::vector<Particle>& particles = sp(cell);
            std::vector<double>& cell_moments = moments(cell);
            cell_moments.assign(particles.size(), 0.0);
            for (size_t i = 0; i < particles.size(); ++i) {
                const Particle& p = particles[i];
                const Vector3R coord = domain.to_cell_coordinates(p.coord);
                const Vector3I origin = checked_floor_index(coord);
                const Stencil3D<2> no = make_linear_stencil_at(coord, origin);
                const Stencil3D<2> sh = make_linear_stencil_at(coord, linear_half_shifted_origin(coord, origin),
                                                               Vector3R::Constant(0.5));
                const Vector3R field = interpolateB(magnetic, no, sh);
                const double strength = field.norm();
                if (!(strength > 0))
                    continue;
                const Vector3R b = field / strength;
                const Vector3R drift = interpolateE(electric, no, sh).cross(field) / (strength * strength);
                Vector3R gyration = p.velocity - p.velocity.dot(b) * b - drift;
                gyration -= gyration.dot(b) * b;
                const double theta = 0.5 * std::abs(qm) * strength * dt;
                cell_moments[i] = 0.5 * mw * gyration.squaredNorm() / strength * theta * theta / (1 + theta * theta);
                const Vector3R dm = -cell_moments[i] * b;
                visit_magnetic_weights(p.coord, domain, [&](int ix, int iy, int iz, int c, double w) {
                    view(ix, iy, iz, c) += w * dm[c];
                });
            }
        }
    }
    for (const Field3d& grid : partial) moment += grid;
}

// dK_p = -dM_p . (B^{n+1} - B^n) at x^{n+1/2} = particle.coord - dt/2 v^{n+1}
// (after the second half-drift, before the wall decisions, particles in the
// order of the deposit); the energy goes into the gyration (v_perp scaled at
// fixed direction). A particle that cannot give it leaves it to the caller.
inline Exchange exchange(ParticlesArray& sp, const Field3d& magnetic, const Field3d& magnetic_change, double dt,
                         const Moments& moments) {
    RECORD_TIMER;
    const Domain& domain = sp.get_domain();
    const double mw = sp.mass() * sp.mpw();
    const GridNodeView<const Field3d> change(magnetic_change);
    double target = 0, given = 0;
#pragma omp parallel for schedule(dynamic, sp.omp_granularity()) reduction(+ : target, given)
    for (int cell = 0; cell < sp.size(); ++cell) {
        std::vector<Particle>& particles = sp(cell);
        const std::vector<double>& cell_moments = moments(cell);
        for (size_t i = 0; i < particles.size(); ++i) {
            if (!(cell_moments[i] > 0))
                continue;
            Particle& p = particles[i];
            const Vector3R middle = p.coord - 0.5 * dt * p.velocity;
            const Vector3R b = field_direction(magnetic, middle, domain);
            const Vector3R dm = -cell_moments[i] * b;
            double work = 0;
            visit_magnetic_weights(middle, domain, [&](int ix, int iy, int iz, int c, double w) {
                work -= w * dm[c] * change(ix, iy, iz, c);
            });
            target += work;
            const Vector3R v_perp = p.velocity - p.velocity.dot(b) * b;
            const double perp2 = v_perp.squaredNorm();
            const double scaled = perp2 + 2 * work / mw;
            if (!(perp2 > 0) || scaled < 0)
                continue;
            p.velocity += (std::sqrt(scaled / perp2) - 1) * v_perp;
            given += work;
        }
    }
    return {target, given};
}

}   // namespace pathMagnetization
