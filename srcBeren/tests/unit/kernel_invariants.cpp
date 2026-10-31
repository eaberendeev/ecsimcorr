#include "algorithms_ecsim.h"
#include "algorithms_ecsim_corr.h"
#include "coupling/current.h"
#include "coupling/density.h"
#include "coupling/gather.h"
#include "coupling/moments.h"
#include "ecsim/response/response.h"
#include "test_support.h"

void yee_gather_reproduces_affine_fields() {
    Domain domain;
    domain.init(Vector3I(6, 5, 4), Vector3R(0.7, 1.1, 0.9));
    Field3d e(domain.size(), 3), b(domain.size(), 3);
    const auto value = [](const Vector3R& r, int d) { return 1.3 + d + 0.2 * r.x() - 0.3 * r.y() + 0.4 * r.z(); };
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                for (int d = 0; d < 3; ++d) {
                    e(i, j, k, d) = value(domain.get_node_position(i, j, k, FieldType::ELECTRIC, d), d);
                    b(i, j, k, d) = value(domain.get_node_position(i, j, k, FieldType::MAGNETIC, d), d);
                }
    for (int sample = 0; sample < 20; ++sample) {
        const Vector3R normalized(1.1 + sample * 0.11, 1.3 + sample * 0.07, 1.2 + sample * 0.03);
        const Vector3R physical(normalized.x() * 0.7, normalized.y() * 1.1, normalized.z() * 0.9);
        ParticleShape<Shape, 2> regular, shifted;
        regular.fill_from_normalized(normalized);
        shifted.fill_from_normalized(normalized, Vector3R(0.5, 0.5, 0.5));
        const auto el = interpolateE_linear(e, normalized), bl = interpolateB_linear(b, normalized);
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
void ecsim_response_matches_midpoint_current_and_work(bool reference) {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(0.7, 1.1, 0.9));
    const double dt = 0.3;
    ParticlesArray sp(check::species("beam", -1, 2, 1, 4), domain);
    sp.add_particle(Particle(Vector3R(1.8, 2.4, 1.5), Vector3R(0.03, -0.02, 0.04)));
    sp.add_particle(Particle(Vector3R(2.4, 3.2, 2.3), Vector3R(-0.01, 0.05, -0.03)));
    Field3d e(domain.size(), 3), b(domain.size(), 3), j0(domain.size(), 3), j1(domain.size(), 3), jp(domain.size(), 3);
    const Vector3R electric(0.02, -0.03, 0.01), magnetic(-0.2, 0.3, 0.1);
    for (int i = 0; i < e.capacity(); ++i) {
        e[i] = electric[i % 3];
        b[i] = magnetic[i % 3];
    }
    j0.setZero();
    j1.setZero();
    jp.setZero();
    const double before = sp.get_kinetic_energy();
    coupling::deposit_point_current(sp, j0);
    algorithmsECSIM::predict_current<ShapeType::Linear>(sp, b, jp, dt);
    Array3D<int> counts(domain.size());
    for (int i = 0; i < sp.size(); ++i) counts(i) = sp.particlesData(i).size();
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
        check::scenario("Linear current: drho/dt + div J = 0 for oblique cell crossings",
                        oblique_trajectories_obey_charge_continuity<ShapeType::Linear>);
        check::scenario("Quadratic current: drho/dt + div J = 0 for oblique cell crossings",
                        oblique_trajectories_obey_charge_continuity<ShapeType::Quadratic>);
        check::scenario("Reference ECSIM gather/Jp/L give midpoint current and particle energy gain",
                        [] { ecsim_response_matches_midpoint_current_and_work(true); });
        check::scenario("Production ECSIM gather/Jp/L give midpoint current and particle energy gain",
                        [] { ecsim_response_matches_midpoint_current_and_work(false); });
    });
}
