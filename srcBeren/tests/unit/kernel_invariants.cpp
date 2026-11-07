#include "coupling/current.h"
#include "coupling/density.h"
#include "coupling/gather.h"
#include "coupling/moments.h"
#include "raw_gather_reference.h"
#include "simulation/ecsim/algorithms_ecsim.h"
#include "simulation/ecsim/response/response.h"
#include "simulation/ecsim_corr/algorithms_ecsim_corr.h"
#include "test_support.h"

void gather_rejects_incompatible_domain_before_particle_loop() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1.0), GHOST_CELLS + 1);
    ParticlesArray species(check::species("beam"), domain);
    Field3d electric(domain.size(), 3), magnetic(domain.size(), 3);
    bool rejected = false;
    try {
        algorithmsECSIM::predict_velocity<ShapeType::Linear>(species, electric, magnetic, 0.1);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "a mismatched gather halo is rejected before OpenMP particle processing");
}

void yee_gather_reproduces_affine_fields() {
    Domain domain;
    domain.init(Vector3I(6, 5, 4), Vector3R(0.7, 1.1, 0.9));
    Field3d e(domain.size(), 3), b(domain.size(), 3);
    const auto value = [](const Vector3R& r, int d) { return 1.3 + d + 0.2 * r.x() - 0.3 * r.y() + 0.4 * r.z(); };
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                for (int d = 0; d < 3; ++d) {
                    e(i, j, k, d) = value(domain.storage_node_position(i, j, k, FieldType::ELECTRIC, d), d);
                    b(i, j, k, d) = value(domain.storage_node_position(i, j, k, FieldType::MAGNETIC, d), d);
                }
    for (int sample = 0; sample < 20; ++sample) {
        const Vector3R normalized(1.1 + sample * 0.11, 1.3 + sample * 0.07, 1.2 + sample * 0.03);
        const Vector3R physical(normalized.x() * 0.7, normalized.y() * 1.1, normalized.z() * 0.9);
        const auto regular = make_point_stencil<ShapeType::Linear>(normalized);
        const auto shifted = make_point_stencil<ShapeType::Linear>(normalized, Vector3R(0.5));
        const auto el = raw_gather_reference::interpolateE_linear(e, normalized),
                   bl = raw_gather_reference::interpolateB_linear(b, normalized);
        check::vector(interpolateE(e, regular, shifted), el, 1e-13, "generic and optimized E interpolation");
        check::vector(interpolateB(b, regular, shifted), bl, 1e-13, "generic and optimized B interpolation");
        for (int d = 0; d < 3; ++d) {
            check::near(el[d], value(physical, d), 1e-13, "Yee E must reproduce an affine field");
            check::near(bl[d], value(physical, d), 1e-13, "Yee B must reproduce an affine field");
        }
    }
}

template <ShapeType ST>
void deposit_preserves_charge_or_particle_number() {
    Domain domain;
    domain.init(Vector3I(6, 5, 4), Vector3R(0.7, 1.1, 0.9));
    for (bool neutral : {false, true}) {
        ParticlesArray sp(check::species("species", neutral ? 0 : -2, 3, 2, 10), domain);
        for (int i = 0; i < 12; ++i)
            sp.add_particle(Particle(Vector3R(0.8 + i * 0.15, 1.4, 1.3), Vector3R(0.03, 0.02, 0.01)));
        Field3d density(domain.size(), 1);
        coupling::deposit_density<ST>(sp, density);
        const Field3d actual = density;
        const double total = actual.data().sum();
        check::near(total, 12 * 0.2 * (neutral ? 1 : -2), 1e-13, "deposition lost total charge/number");
        coupling::deposit_density_reference<ST>(sp, density);
        check::near((actual - density).norm(), 0, 1e-13, "blocked density differs from reference");
    }
}

// Partition of unity gives the integral current q*mpw*sum(v). Depositing
// twice must accumulate; depositing density/second moment must overwrite.
void common_deposits_have_explicit_accumulation_contracts() {
    Domain domain;
    domain.init(Vector3I(6, 5, 4), Vector3R(0.7, 1.1, 0.9));
    ParticlesArray species(check::species("beam", -2, 3, 2, 10), domain);
    const Vector3R v1(1, -2, 3), v2(-2, 1, 0.5);
    species.add_particle(Particle(Vector3R(1.8, 2.4, 1.5), v1));
    species.add_particle(Particle(Vector3R(2.4, 3.2, 2.3), v2));
    Field3d current(domain.size(), 3);
    current.setZero();
    const Vector3R expected_current = -0.4 * (v1 + v2);
    for (int repeat = 1; repeat <= 2; ++repeat) {
        coupling::deposit_point_current(species, current);
        Vector3R integral(0, 0, 0);
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k)
                    for (int d = 0; d < 3; ++d) integral[d] += current(i, j, k, d);
        check::vector(integral, double(repeat) * expected_current, 1e-13, "point current accumulates q*mpw*v");
    }
    Field3d moment(domain.size(), 1);
    moment.data().setConstant(17);
    const auto vx = [](const Vector3R&, const Vector3R& v, const Vector3R&) { return v.x(); };
    const auto vy = [](const Vector3R&, const Vector3R& v, const Vector3R&) { return v.y(); };
    coupling::deposit_second_moment(species, moment, vx, vy);
    check::near(moment.data().sum(), 0.6 * (v1.x() * v1.y() + v2.x() * v2.y()), 1e-13,
                "mixed second moment overwrites and retains the sign of vx*vy");
    moment.data().setConstant(17);
    coupling::deposit_density<ShapeType::Linear>(species, moment);
    check::near(moment.data().sum(), -0.8, 1e-13, "density overwrites with charge-weighted shape");
}

void kinetic_moments_are_independent_of_threads() {
    Domain domain;
    domain.init(Vector3I(8, 6, 4), Vector3R(1, 1, 1));
    ParticlesArray sp(check::species("beam", -1, 2, 1, 1), domain);
    long double expected[3] = {};
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 6; ++j)
            for (int k = 0; k < 4; ++k) {
                const Vector3R velocity(0.01 * (1 + i), 0.03 * (1 + j), 0.07 * (1 + k));
                sp.add_particle(Particle(Vector3R(i + 0.25, j + 0.25, k + 0.25), velocity));
                for (int d = 0; d < 3; ++d) expected[d] += static_cast<long double>(velocity[d]) * velocity[d];
            }
    const int original_threads = omp_get_max_threads();
    omp_set_num_threads(1);
    const ParticleMoments baseline = sp.diagnostic_moments();
    for (int threads : {2, 4}) {
        omp_set_num_threads(threads);
        const ParticleMoments actual = sp.diagnostic_moments();
        check::equal(actual.count, baseline.count, "moment count depends on thread count");
        check::vector(actual.kinetic, baseline.kinetic, 0, "moment sums depend on thread count");
    }
    omp_set_num_threads(original_threads);
    check::equal(baseline.count, 8 * 6 * 4, "moments lost particles");
    for (int d = 0; d < 3; ++d)
        check::near(baseline.kinetic[d], static_cast<double>(expected[d]), 1e-13,
                    "kinetic moment differs from the independent particle sum");
}

void cylindrical_second_moment_has_correct_energy_trace() {
    for (const Vector3I cells : {Vector3I(6, 6, 6), Vector3I(6, 4, 6)}) {
        Domain domain;
        domain.init(cells, Vector3R(1, 1, 1));
        const Vector3R center(0.5 * cells.x(), 0.5 * cells.y(), 2.0);
        for (const Vector3R offset : {Vector3R(0, 0, 0), Vector3R(0.3, 0.4, 0)}) {
            ParticlesArray sp(check::species("beam", -1, 2, 1, 1), domain);
            sp.add_particle(Particle(center + offset, Vector3R(1, 2, 3)));
            Field3d rr(domain.size(), 1), pp(domain.size(), 1), zz(domain.size(), 1);
            coupling::deposit_second_moment(sp, rr, coupling::radial_velocity, coupling::radial_velocity);
            coupling::deposit_second_moment(sp, pp, coupling::phi_velocity, coupling::phi_velocity);
            coupling::deposit_second_moment(sp, zz, coupling::z_velocity, coupling::z_velocity);
            const bool axis = offset.x() == 0;
            const double vr = axis ? 1 : 2.2, vp = axis ? 2 : 0.4;
            check::near(rr.data().sum(), 2 * vr * vr, 1e-12,
                        "radial pressure uses the physical center and finite axis basis");
            check::near(pp.data().sum(), 2 * vp * vp, 1e-12, "azimuthal pressure");
            check::near(rr.data().sum() + pp.data().sum() + zz.data().sum(), 2 * sp.get_kinetic_energy(), 1e-12,
                        "pressure trace must equal twice kinetic energy");
        }
    }
}

// Raw placement depends on the actual Grid halo. Check every stored node for
// a complete footprint and for clipping at either edge, with nonuniform data
// so an incorrect origin or transposed axis cannot preserve the expectation.
void current_buffer_flush_uses_actual_grid_footprint() {
    for (int ghost : {0, 1, 2}) {
        const Grid grid(Vector3R(0.7, 1.1, 0.9), Vector3I(4, 5, 6), ghost);
        const Vector3I sizes = grid.size();
        const std::array<Vector3I, 5> raw_origins{{Vector3I(0, 1, 2), Vector3I(-1, 1, 2), sizes - Vector3I(3, 2, 1),
                                                   Vector3I(-5, 0, 0), Vector3I(sizes.x() + 1, 0, 0)}};
        CurrentBuffer<4> buffer;
        for (int x = 0; x < 4; ++x)
            for (int y = 0; y < 4; ++y)
                for (int z = 0; z < 4; ++z)
                    for (int d = 0; d < 3; ++d) buffer(x, y, z, d) = 1 + x + 4 * y + 16 * z + 64 * d;
        for (const auto& raw_origin : raw_origins)
            for (bool checked : {false, true}) {
                Field3d field(sizes, 3);
                field.data().setConstant(0.125);
                const Vector3I logical_origin = raw_origin - Vector3I(ghost);
                if (checked)
                    flush_current_buffer_with_check(field, buffer,
                                                    to_storage_index(logical_origin, grid.ghost_cells()));
                else
                    flush_current_buffer(field, buffer, to_storage_index(logical_origin, grid.ghost_cells()));
                for (int i = 0; i < sizes.x(); ++i)
                    for (int j = 0; j < sizes.y(); ++j)
                        for (int k = 0; k < sizes.z(); ++k) {
                            const Vector3I local = Vector3I(i, j, k) - raw_origin;
                            const bool inside = local.x() >= 0 && local.x() < 4 && local.y() >= 0 && local.y() < 4 &&
                                                local.z() >= 0 && local.z() < 4;
                            for (int d = 0; d < 3; ++d) {
                                const double expected =
                                    0.125 + (inside ? buffer(local.x(), local.y(), local.z(), d) : 0);
                                check::near(field(i, j, k, d), expected, 0,
                                            "current flush matches raw footprint, ghost=" + std::to_string(ghost));
                            }
                        }
            }
    }
}

// Include populated lower halo bins (raw bin zero), upper halo bins, and an
// interior crossing. All nonzero endpoint weights fit in the allocated halo;
// envelope padding can extend beyond it. This is a deposit test without wall
// handling, so the mathematically prescribed trajectories can lie in the halo.
template <ShapeType ST>
void trajectory_current_uses_actual_halo() {
    for (int ghost : {0, 1, 2}) {
        Domain domain;
        const Vector3R cell_size(0.7, 1.1, 0.9);
        domain.init(Vector3I(8, 7, 6), cell_size, ghost);
        ParticlesArray optimized(check::species("beam", -1, 1, 1, 1), domain);
        ParticlesArray reference(check::species("beam", -1, 1, 1, 1), domain);
        const std::array<Vector3R, 3> starts{{Vector3R(0.75 - ghost, 0.8 - ghost, 0.85 - ghost),
                                              Vector3R(8 + ghost - 0.75, 7 + ghost - 0.8, 6 + ghost - 0.85),
                                              Vector3R(2.9, 3.1, 2.8)}};
        const std::array<Vector3R, 3> velocities{
            {Vector3R(0.1, 0.15, 0.05), Vector3R(0.1, 0.15, 0.05), Vector3R(0.8, -0.7, 0.4)}};
        Vector3R expected_current(0.0);
        for (size_t n = 0; n < starts.size(); ++n) {
            Vector3R position, velocity;
            for (int d = 0; d < 3; ++d) {
                position[d] = starts[n][d] * cell_size[d];
                velocity[d] = velocities[n][d] * cell_size[d];
            }
            const Particle particle(position, velocity);
            optimized.add_particle(particle);
            reference.add_particle(particle);
            expected_current -= velocity;
        }
        check::equal(optimized.cells()(0, 0, 0).size(), size_t(1), "fixture populates the lower storage bin");
        Field3d density(domain.size(), 1);
        coupling::deposit_density<ST>(optimized, density);
        const Field3d before = density;
        optimized.currentOnGrid.setZero();
        reference.currentOnGrid.setZero();
        coupling::drift_and_deposit_current<ST>(optimized, 1, optimized.currentOnGrid);
        coupling::drift_and_deposit_current_reference<ST>(reference, 1, reference.currentOnGrid);
        optimized.update_cells(domain);
        coupling::deposit_density<ST>(optimized, density);
        Vector3R total_current(0.0);
        const std::string context = "trajectory deposit, ghost=" + std::to_string(ghost);
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k) {
                    const auto& current = optimized.currentOnGrid;
                    const double div = (current(i, j, k, 0) - (i > 0 ? current(i - 1, j, k, 0) : 0)) / cell_size.x() +
                                       (current(i, j, k, 1) - (j > 0 ? current(i, j - 1, k, 1) : 0)) / cell_size.y() +
                                       (current(i, j, k, 2) - (k > 0 ? current(i, j, k - 1, 2) : 0)) / cell_size.z();
                    check::near(density(i, j, k, 0) - before(i, j, k, 0) + div, 0, 2e-13,
                                context + " charge continuity on every stored node");
                    for (int d = 0; d < 3; ++d) {
                        check::near(current(i, j, k, d), reference.currentOnGrid(i, j, k, d), 2e-13,
                                    context + " optimized/reference agreement");
                        total_current[d] += current(i, j, k, d);
                    }
                }
        check::vector(total_current, expected_current, 2e-13, context + " total current equals q*mpw*displacement/dt");
    }
}

void second_moment_uses_actual_halo() {
    for (int ghost : {1, 2}) {
        Domain domain;
        const Vector3R cell_size(0.7, 1.1, 0.9);
        domain.init(Vector3I(5, 6, 7), cell_size, ghost);
        ParticlesArray species(check::species("moment", -1, 3, 2, 4), domain);
        const std::array<Vector3R, 2> positions{
            {Vector3R(0.25 - ghost, 0.375 - ghost, 0.625 - ghost), Vector3R(2.3, 3.4, 4.2)}};
        double expected_total = 0;
        Vector3R expected_first(0.0);
        for (size_t n = 0; n < positions.size(); ++n) {
            Vector3R position;
            for (int d = 0; d < 3; ++d) position[d] = positions[n][d] * cell_size[d];
            const Vector3R velocity(1.5 + n, -2, 0.75);
            species.add_particle(Particle(position, velocity));
            const double moment = 3 * 0.5 * velocity.x() * velocity.x();
            expected_total += moment;
            expected_first += moment * position;
        }
        Field3d moment(domain.size(), 1);
        moment.data().setConstant(19);
        const auto vx = [](const Vector3R&, const Vector3R& velocity, const Vector3R&) { return velocity.x(); };
        coupling::deposit_second_moment(species, moment, vx, vx);
        double actual_total = 0;
        Vector3R actual_first(0.0);
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k) {
                    const double value = moment(i, j, k, 0);
                    check::require(std::isfinite(value), "finite second-moment deposit with negative lower halo point");
                    actual_total += value;
                    actual_first += value * Vector3R((i - ghost) * cell_size.x(), (j - ghost) * cell_size.y(),
                                                     (k - ghost) * cell_size.z());
                }
        check::near(actual_total, expected_total, 1e-13, "second-moment total for actual grid halo");
        check::vector(actual_first, expected_first, 1e-13, "second-moment spatial first moment for actual grid halo");
    }
}

void second_moment_skips_zero_padding_at_upper_storage_node() {
    Domain domain;
    domain.init(Vector3I(2, 2, 2), Vector3R(0.5, 1, 2), 2);
    ParticlesArray species(check::species("moment", -1, 3, 1, 1), domain);
    species.add_particle(Particle(Vector3R(2, 4, 8), Vector3R(2, 0, 0)));
    Field3d moment(domain.size(), 1);
    const auto vx = [](const Vector3R&, const Vector3R& velocity, const Vector3R&) { return velocity.x(); };
    coupling::deposit_second_moment(species, moment, vx, vx);
    check::near(moment.data().sum(), 12, 1e-14, "upper-node moment retains its full integral");
    check::near(moment(domain.size().x() - 1, domain.size().y() - 1, domain.size().z() - 1, 0), 12, 1e-14,
                "zero support padding beyond allocated storage is not accessed");
}

// Endpoints share four logical nodes on every axis, including when the path
// crosses a cell boundary. Compare every component, including zero padding and
// translated buffers; independent charge continuity is checked below.
template <ShapeType ST>
void esirkepov_stencil_kernels_agree() {
    constexpr auto weight = ShapeKernel<ST>::weight;
    const std::array<std::pair<Vector3R, Vector3R>, 5> segments{{
        {Vector3R(2.9, 3.1, 3.8), Vector3R(3.9, 2.1, 4.5)},
        {Vector3R(3.1, 3.1, 3.8), Vector3R(2.1, 2.1, 4.5)},
        {Vector3R(-0.5, -0.1, 0.5), Vector3R(0.5, -1.1, 1.2)},
        {Vector3R(0.5, 1.0, -0.5), Vector3R(0.5, 1.0, -0.5)},
        {Vector3R(0.499999, 0.5, 0.500001), Vector3R(-0.500001, 1.5, -0.499999)},
    }};
    for (const auto& [from, to] : segments) {
        Vector3I origin;
        for (int d = 0; d < 3; ++d)
            origin[d] =
                static_cast<int>(std::floor(std::min(from[d], to[d]) + (ST == ShapeType::Quadratic ? 0.5 : 0.0))) - 1;
        const auto start = make_stencil_at<weight, 4>(from, origin);
        const auto end = make_stencil_at<weight, 4>(to, origin);
        CurrentBuffer<4> reference, optimized;
        optimized.zero();
        decompose_esirkepov_current_reference(start, end, 0.7, -1.1, 0.9, reference);
        decompose_esirkepov_current(start, end, 0.7, -1.1, 0.9, optimized);
        CurrentBuffer<7> translated;
        translated.data_.fill(0.25);
        decompose_esirkepov_current_offsetted(start, end, 0.7, -1.1, 0.9, translated, 1, 2, 3);
        for (int x = 0; x < 7; ++x)
            for (int y = 0; y < 7; ++y)
                for (int z = 0; z < 7; ++z)
                    for (int d = 0; d < 3; ++d) {
                        const bool inside = x >= 1 && x < 5 && y >= 2 && y < 6 && z >= 3 && z < 7;
                        const double expected = inside ? reference(x - 1, y - 2, z - 3, d) : 0;
                        check::require(std::isfinite(translated(x, y, z, d)), "finite translated Esirkepov current");
                        check::near(translated(x, y, z, d), 0.25 + expected, 1e-14,
                                    "translated Esirkepov kernel accumulates into the correct nodes");
                    }
        for (size_t i = 0; i < reference.data_.size(); ++i) {
            check::require(std::isfinite(reference.data_[i]) && std::isfinite(optimized.data_[i]),
                           "finite Esirkepov stencil current");
            check::near(optimized.data_[i], reference.data_[i], 1e-14,
                        "optimized Esirkepov stencil kernel agrees with reference");
        }
    }
}

template <ShapeType ST>
void oblique_trajectories_obey_charge_continuity() {
    Domain domain;
    const Vector3R cell_size(0.7, 1.1, 0.9);
    domain.init(Vector3I(8, 8, 8), cell_size);
    for (double displacement : {0.8, 1.0})
        for (int direction : {-1, 1})
            for (int impl = 0; impl < 3; ++impl) {
                ParticlesArray sp(check::species("beam"), domain);
                sp.add_particle(Particle(
                    Vector3R((direction > 0 ? 2.9 : 3.1) * cell_size.x(), 3.1 * cell_size.y(), 3.8 * cell_size.z()),
                    Vector3R(direction * displacement * cell_size.x(), -displacement * cell_size.y(),
                             0.7 * cell_size.z())));
                Field3d density(domain.size(), 1);
                coupling::deposit_density<ST>(sp, density);
                const Field3d before = density;
                sp.currentOnGrid.setZero();
                if (impl == 0)
                    coupling::drift_and_deposit_current<ST>(sp, 1.0, sp.currentOnGrid);
                else if (impl == 1)
                    coupling::drift_and_deposit_current_reference<ST>(sp, 1.0, sp.currentOnGrid);
                else {
                    Field3d zero(domain.size(), 3);
                    zero.setZero();
                    double work = 0;
                    fused_push_and_deposit<ST>(sp, zero, zero, 2.0, 1.0, sp.currentOnGrid, work);
                }
                sp.update_cells(domain);
                coupling::deposit_density<ST>(sp, density);
                double error = 0;
                for (int i = 1; i < domain.size().x() - 1; ++i)
                    for (int j = 1; j < domain.size().y() - 1; ++j)
                        for (int k = 1; k < domain.size().z() - 1; ++k) {
                            const auto& current = sp.currentOnGrid;
                            const double div = (current(i, j, k, 0) - current(i - 1, j, k, 0)) / cell_size.x() +
                                               (current(i, j, k, 1) - current(i, j - 1, k, 1)) / cell_size.y() +
                                               (current(i, j, k, 2) - current(i, j, k - 1, 2)) / cell_size.z();
                            const double residual = density(i, j, k, 0) - before(i, j, k, 0) + div;
                            check::require(std::isfinite(residual), "nonfinite moving-charge continuity residual");
                            error = std::max(error, std::abs(residual));
                        }
                check::near(error, 0, 1e-13, "charge continuity for an oblique trajectory up to one cell per axis");
            }
}
void ecsim_response_matches_midpoint_current_and_work(bool reference, bool varying_fields = false) {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(0.7, 1.1, 0.9));
    const double dt = 0.3;
    ParticlesArray sp(check::species("beam", -1, 2, 1, 4), domain);
    sp.add_particle(Particle(Vector3R(1.8, 2.4, 1.5), Vector3R(0.03, -0.02, 0.04)));
    sp.add_particle(Particle(Vector3R(2.4, 3.2, 2.3), Vector3R(-0.01, 0.05, -0.03)));
    Field3d e(domain.size(), 3), b(domain.size(), 3), j0(domain.size(), 3), j1(domain.size(), 3), jp(domain.size(), 3);
    const Vector3R electric(0.02, -0.03, 0.01), magnetic(-0.2, 0.3, 0.1);
    for (int i = 0; i < e.capacity(); ++i) {
        e[i] = electric[i % 3] + (varying_fields ? 0.007 * std::sin(0.37 * i) : 0.0);
        b[i] = magnetic[i % 3] + (varying_fields ? 0.03 * std::cos(0.23 * i) : 0.0);
    }
    j0.setZero();
    j1.setZero();
    jp.setZero();
    const double before = sp.get_kinetic_energy();
    coupling::deposit_point_current(sp, j0);
    algorithmsECSIM::predict_current<ShapeType::Linear>(sp, b, jp, dt);
    Array3D<int> counts(domain.size());
    for (int i = 0; i < sp.size(); ++i) counts(i) = sp(i).size();
    ecsim_response::Workspace workspace;
    Operator response(3 * domain.total_size(), 3 * domain.total_size());
    if (reference) {
        BlockMatrix reference_blocks;
        reference_blocks.resize(domain.total_size());
        reference_blocks.prepare(counts);
        ecsim_response::accumulate_species_reference<ShapeType::Linear>(sp, b, domain, dt, reference_blocks);
        ecsim_response::assemble_reference(response, domain, reference_blocks);
    } else {
        ecsim_response::RowBlocks row_blocks;
        ecsim_response::clear_row_blocks(row_blocks);
        ecsim_response::accumulate_species<ShapeType::Linear>(sp, b, domain, dt, row_blocks);
        ecsim_response::assemble_matrix(response, row_blocks, workspace);
    }
    algorithmsECSIM::predict_velocity<ShapeType::Linear>(sp, e, b, dt);
    coupling::deposit_point_current(sp, j1);
    const Eigen::VectorXd mid = 0.5 * (j0.data() + j1.data());
    const Eigen::VectorXd predicted = jp.data() + (2.0 / dt) * (response * e.data());
    check::near((mid - predicted).norm(), 0, 1e-13,
                "predicted current/mass matrix disagrees with actual midpoint current");
    check::near(sp.get_kinetic_energy() - before, dt * e.data().dot(mid), 1e-13, "gather/scatter work is not adjoint");
}

int main() {
    return check::run([] {
        check::scenario("Gather rejects mismatched compile-time halo",
                        gather_rejects_incompatible_domain_before_particle_loop);
        check::scenario("Second moment at the last storage node",
                        second_moment_skips_zero_padding_at_upper_storage_node);
        check::scenario("Common deposits accumulate current and overwrite density and second moments",
                        common_deposits_have_explicit_accumulation_contracts);
        check::scenario("Yee gather reproduces a known affine E/B field", yee_gather_reproduces_affine_fields);
        check::scenario("Linear deposit preserves total charge and neutral particle number",
                        deposit_preserves_charge_or_particle_number<ShapeType::Linear>);
        check::scenario("Quadratic deposit preserves total charge and neutral particle number",
                        deposit_preserves_charge_or_particle_number<ShapeType::Quadratic>);
        check::scenario("Particle kinetic moments agree with a serial sum and do not depend on threads",
                        kinetic_moments_are_independent_of_threads);
        check::scenario("Cylindrical second moment has trace 2*K (bulk flow is included)",
                        cylindrical_second_moment_has_correct_energy_trace);
        check::scenario("Current-buffer flush places and clips its full footprint for halos 0/1/2",
                        current_buffer_flush_uses_actual_grid_footprint);
        check::scenario("Linear current preserves continuity and q*v across populated halo bins 0/1/2",
                        trajectory_current_uses_actual_halo<ShapeType::Linear>);
        check::scenario("Quadratic current preserves continuity and q*v across populated halo bins 0/1/2",
                        trajectory_current_uses_actual_halo<ShapeType::Quadratic>);
        check::scenario("Second moment preserves total and spatial first moment in halos 1/2",
                        second_moment_uses_actual_halo);
        check::scenario("Linear Esirkepov kernels agree for shared and translated four-node stencils",
                        esirkepov_stencil_kernels_agree<ShapeType::Linear>);
        check::scenario("Quadratic Esirkepov kernels agree for shared and translated four-node stencils",
                        esirkepov_stencil_kernels_agree<ShapeType::Quadratic>);
        check::scenario("Linear current: drho/dt + div J = 0 for oblique cell crossings",
                        oblique_trajectories_obey_charge_continuity<ShapeType::Linear>);
        check::scenario("Quadratic current: drho/dt + div J = 0 for oblique cell crossings",
                        oblique_trajectories_obey_charge_continuity<ShapeType::Quadratic>);
        check::scenario("Reference ECSIM gather/Jp/L give midpoint current and particle energy gain",
                        [] { ecsim_response_matches_midpoint_current_and_work(true); });
        check::scenario("Production ECSIM gather/Jp/L give midpoint current and particle energy gain",
                        [] { ecsim_response_matches_midpoint_current_and_work(false); });
        check::scenario("Spatially varying fields: reference ECSIM response matches midpoint current and work",
                        [] { ecsim_response_matches_midpoint_current_and_work(true, true); });
        check::scenario("Spatially varying fields: production ECSIM response matches midpoint current and work",
                        [] { ecsim_response_matches_midpoint_current_and_work(false, true); });
    });
}
