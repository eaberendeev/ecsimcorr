#include "ParticlesDiagnostic.h"
#include "algorithms_ecsim.h"
#include "algorithms_ecsim_corr.h"
#include "env_options.h"
#include "interpolation.h"
#include "test_support.h"

void environment_options() {
    const char* name = "BEREN_TEST_ENV_PARSER";
    check::equal(envOptions::detail::parseIntegerEnv<int>(name, "+42"), 42, "integer plus sign");
    check::equal(envOptions::detail::parseIntegerEnv<int>(name, "-42"), -42, "integer minus sign");
    check::equal(envOptions::detail::parseIntegerEnv<int64_t>(name, "9223372036854775807"),
                 std::numeric_limits<int64_t>::max(), "int64 maximum");
    for (const char* value : {"", " ", "+", "+-1", "++1", "1x", " 1", "1 ", "1\t", "2147483648", "-2147483649"}) {
        bool rejected = false;
        try {
            envOptions::detail::parseIntegerEnv<int>(name, value);
        } catch (const std::invalid_argument& error) {
            rejected = std::string(error.what()).find(name) != std::string::npos;
        }
        check::require(rejected, "invalid integer must identify its environment variable");
    }
    const char* previous = std::getenv(name);
    const bool was_set = previous != nullptr;
    const std::string saved = was_set ? previous : "";
    for (const char* value : {"true", "YES", "1", "false", "No", "0"}) {
        setenv(name, value, 1);
        const bool expected = std::string(value) == "true" || std::string(value) == "YES" || std::string(value) == "1";
        check::equal(envOptions::getenvParsed<bool>(name, !expected), expected, "boolean parsing");
    }
    for (const char* value : {"", " true", "false ", "1\t", "maybe"}) {
        setenv(name, value, 1);
        bool rejected = false;
        try {
            envOptions::getenvParsed<bool>(name, false);
        } catch (const std::invalid_argument& error) {
            rejected = std::string(error.what()).find(name) != std::string::npos;
        }
        check::require(rejected, "invalid boolean must identify its environment variable");
    }
    unsetenv(name);
    check::require(envOptions::getenvParsed<bool>(name, true), "unset boolean uses default");
    check::equal(envOptions::getenvParsed<int>(name, 17), 17, "unset integer uses default");
    if (was_set)
        setenv(name, saved.c_str(), 1);
}

void pusher() {
    const Vector3R initial(0.2, -0.1, 0.4), electric(0.03, -0.02, 0.01), magnetic(0, 0, 0.7);
    Particle p(Vector3R(0, 0, 0), initial);
    borisPusher::update_vEB(p, -0.5, electric, Vector3R(0, 0, 0), 0.4);
    check::vector(p.velocity, initial - 0.2 * electric, 1e-15, "Boris electric impulse");
    p.velocity = Vector3R(1, 0, 0);
    borisPusher::update_vEB(p, 1, Vector3R(0, 0, 0), magnetic, 0.1);
    check::require(p.velocity.y() < 0, "Lorentz rotation has wrong sign");
    p.velocity = initial;
    for (int step = 0; step < 2000; ++step) borisPusher::update_vEB(p, -0.5, Vector3R(0, 0, 0), magnetic, 0.1);
    check::near(p.velocity.squared(), initial.squared(), 1e-13, "magnetic rotation changed kinetic energy");
    p.velocity = initial;
    borisPusher::update_vEB(p, -0.5, electric, magnetic, 0.4);
    borisPusher::update_vEB(p, -0.5, electric, magnetic, -0.4);
    check::vector(p.velocity, initial, 1e-14, "Boris reversibility");
}

void interpolation() {
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
void density() {
    Domain domain;
    domain.init(Vector3I(6, 5, 4), Vector3R(0.7, 1.1, 0.9));
    for (bool neutral : {false, true}) {
        ParticlesArray sp(check::species("species", neutral ? 0 : -2, 3, 2, 10), domain);
        for (int i = 0; i < 12; ++i)
            sp.add_particle(Particle(Vector3R(0.8 + i * 0.15, 1.4, 1.3), Vector3R(0.03, 0.02, 0.01)));
        sp.density_on_grid_update<ST>();
        const Field3d actual = sp.densityOnGrid;
        const double total = actual.data().sum();
        check::near(total, 12 * 0.2 * (neutral ? 1 : -2), 1e-13, "deposition lost total charge/number");
        sp.density_on_grid_update_reference<ST>();
        check::near((actual - sp.densityOnGrid).norm(), 0, 1e-13, "blocked density differs from reference");
    }
}

void reproducible_moments() {
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

void pressure() {
    for (const Vector3I cells : {Vector3I(6, 6, 6), Vector3I(6, 4, 6)}) {
        Domain domain;
        domain.init(cells, Vector3R(1, 1, 1));
        const Vector3R center(0.5 * cells.x(), 0.5 * cells.y(), 2.0);
        for (const Vector3R offset : {Vector3R(0, 0, 0), Vector3R(0.3, 0.4, 0)}) {
            ParticlesArray sp(check::species("beam", -1, 2, 1, 1), domain);
            sp.add_particle(Particle(center + offset, Vector3R(1, 2, 3)));
            Field3d rr(domain.size(), 1), pp(domain.size(), 1), zz(domain.size(), 1);
            sp.calculate_pressure_component(rr, radial_velocity, radial_velocity);
            sp.calculate_pressure_component(pp, phi_velocity, phi_velocity);
            sp.calculate_pressure_component(zz, z_velocity, z_velocity);
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
void moving_charge() {
    Domain domain;
    const Vector3R cell_size(0.7, 1.1, 0.9);
    domain.init(Vector3I(8, 8, 8), cell_size);
    for (int direction : {-1, 1})
        for (int impl = 0; impl < 3; ++impl) {
            ParticlesArray sp(check::species("beam"), domain);
            sp.add_particle(Particle(
                Vector3R((direction > 0 ? 2.9 : 3.1) * cell_size.x(), 3.1 * cell_size.y(), 3.8 * cell_size.z()),
                Vector3R(direction * 0.8 * cell_size.x(), -0.8 * cell_size.y(), 0.7 * cell_size.z())));
            sp.density_on_grid_update<ST>();
            const Field3d before = sp.densityOnGrid;
            sp.currentOnGrid.setZero();
            if (impl == 0)
                move_and_calc_current<ST>(sp, 1.0, sp.currentOnGrid);
            else if (impl == 1)
                move_and_calc_current_reference<ST>(sp, 1.0, sp.currentOnGrid);
            else {
                Field3d zero(domain.size(), 3);
                zero.setZero();
                double work = 0;
                fused_push_and_deposit<ST>(sp, zero, zero, 2.0, 1.0, sp.currentOnGrid, work);
            }
            sp.update_cells(domain);
            sp.density_on_grid_update<ST>();
            double error = 0;
            for (int i = 1; i < domain.size().x() - 1; ++i)
                for (int j = 1; j < domain.size().y() - 1; ++j)
                    for (int k = 1; k < domain.size().z() - 1; ++k) {
                        const auto& current = sp.currentOnGrid;
                        const double div = (current(i, j, k, 0) - current(i - 1, j, k, 0)) / cell_size.x() +
                                           (current(i, j, k, 1) - current(i, j - 1, k, 1)) / cell_size.y() +
                                           (current(i, j, k, 2) - current(i, j, k - 1, 2)) / cell_size.z();
                        error = std::max(error, std::abs(sp.densityOnGrid(i, j, k, 0) - before(i, j, k, 0) + div));
                    }
            check::near(error, 0, 1e-13, "charge continuity for an oblique 0.8-cell trajectory");
        }
}
void response_work() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(0.7, 1.1, 0.9));
    BoundaryConditionHandler bc;
    Mesh mesh;
    const double dt = 0.3;
    mesh.init(domain, dt, bc);
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
    algorithmsECSIM::calculate_current(sp, j0);
    algorithmsECSIM::predict_current<ShapeType::Linear>(sp, b, jp, dt);
    Array3D<int> counts(domain.size());
    for (int i = 0; i < sp.size(); ++i) counts(i) = sp.particlesData(i).size();
    mesh.LmatX2.prepare(counts);
    sp.fill_matrixL2_Reference<ShapeType::Linear>(mesh, b, domain, dt);
    Operator response(3 * domain.total_size(), 3 * domain.total_size());
    mesh.stencil_Lmat2_Reference(response, domain);
    algorithmsECSIM::predict_velocity<ShapeType::Linear>(sp, e, b, dt);
    algorithmsECSIM::calculate_current(sp, j1);
    const Eigen::VectorXd mid = 0.5 * (j0.data() + j1.data());
    const Eigen::VectorXd predicted = jp.data() + (2.0 / dt) * (response * e.data());
    check::near((mid - predicted).norm(), 0, 1e-13,
                "predicted current/mass matrix disagrees with actual midpoint current");
    check::near(sp.get_kinetic_energy() - before, dt * e.data().dot(mid), 1e-13, "gather/scatter work is not adjoint");
}

int main() {
    return check::run([] {
        environment_options();
        pusher();
        interpolation();
        density<ShapeType::Linear>();
        density<ShapeType::Quadratic>();
        reproducible_moments();
        pressure();
        moving_charge<ShapeType::Linear>();
        moving_charge<ShapeType::Quadratic>();
        response_work();
    });
}
