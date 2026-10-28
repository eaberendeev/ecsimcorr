#include "simulation_ecsim.h"
#include "simulation_ecsim_corr.h"
#include "test_support.h"

namespace {

nlohmann::json settings() {
    return {{"Scheme", "ecsim"},
            {"Dt", 0.2},
            {"NumCellsX", 3},
            {"NumCellsY", 3},
            {"NumCellsZ", 3},
            {"Dx", 1},
            {"Dy", 1},
            {"Dz", 1},
            {"StartFromTime", 0},
            {"StartTimeStep", 0},
            {"LastTimestep", 1},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collider", "None"},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

nlohmann::json settings_corr() {
    auto config = settings();
    config["Scheme"] = "ecsim_corr";
    return config;
}

void fill_constant(Field3d& field, const Vector3R& value) {
    const auto sizes = field.sizes();
    for (int i = 0; i < sizes[0]; ++i)
        for (int j = 0; j < sizes[1]; ++j)
            for (int k = 0; k < sizes[2]; ++k)
                for (int d = 0; d < 3; ++d) field(i, j, k, d) = value[d];
}

void test_constant_fields_with_zero_lmat2() {
    char program[] = "field_solve";
    char* argv[] = {program};
    SimulationEcsim simulation(settings(), {{"particles", nlohmann::json::array()}}, 1, argv);
    simulation.init();

    const Vector3R e_value(0.7, -0.2, 0.4);
    const Vector3R b_value(0.3, -0.1, 0.5);
    const Vector3R j_value(0.4, 0.1, -0.3);
    fill_constant(simulation.fieldE, e_value);
    fill_constant(simulation.fieldB, b_value);
    fill_constant(simulation.fieldJp, j_value);
    simulation.fieldEp.setZero();

    // No particles were added, so the correction matrix has no entries.
    check::near(simulation.mesh.Lmat2.norm(), 0.0, 0.0, "empty plasma must leave Lmat2 zero");

    // A constant B has zero discrete curl; build that term explicitly instead
    // of assuming it vanishes so the production 4-argument RHS is exercised.
    Field3d curl_b(simulation.domain.size(), 3);
    curl_b = simulation.mesh.curlB * simulation.fieldB;
    check::near(curl_b.norm(), 0.0, 1e-14, "constant B must have zero discrete curl");

    simulation.predict_electric_field(simulation.fieldEp, simulation.fieldE, simulation.fieldB, simulation.fieldJp);

    const auto sizes = simulation.fieldEp.sizes();
    for (int i = 0; i < sizes[0]; ++i)
        for (int j = 0; j < sizes[1]; ++j)
            for (int k = 0; k < sizes[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const double expected =
                        e_value[d] - 0.5 * simulation.dt_ * j_value[d] + 0.5 * simulation.dt_ * curl_b(i, j, k, d);
                    check::near(simulation.fieldEp(i, j, k, d), expected, 1e-11,
                                "zero-Lmat2 predictor must give Ep = E + 0.5*dt*(curlB - J)");
                }

    const Field3d rhs = simulation.fieldE - 0.5 * simulation.dt_ * simulation.fieldJp + 0.5 * simulation.dt_ * curl_b;
    const double residual = (simulation.mesh.IMmat * simulation.fieldEp.data() - rhs.data()).norm() / rhs.data().norm();
    check::near(residual, 0.0, 1e-11, "zero-Lmat2 solve must satisfy IMmat*Ep = E + 0.5*dt*(curlB - J)");
}

void test_synthetic_diagonal_lmat2() {
    char program[] = "field_solve";
    char* argv[] = {program};
    SimulationEcsim simulation(settings(), {{"particles", nlohmann::json::array()}}, 1, argv);
    simulation.init();

    const Vector3R e_value(0.7, -0.2, 0.4);
    const Vector3R b_value(0.3, -0.1, 0.5);
    const Vector3R j_value(0.4, 0.1, -0.3);
    fill_constant(simulation.fieldE, e_value);
    fill_constant(simulation.fieldB, b_value);
    fill_constant(simulation.fieldJp, j_value);
    simulation.fieldEp.setZero();

    // Synthetic algebraic check of the A = IMmat + Lmat2 composition: a constant
    // Ep makes IMmat*Ep = Ep (curlE annihilates constants) and Lmat2*Ep =
    // alpha*Ep. This deliberately does not test the assembled physical L, which
    // is covered by lmat2_assembly.
    Operator identity(simulation.mesh.IMmat.rows(), simulation.mesh.IMmat.cols());
    identity.setIdentity();
    const double alpha = 0.25;
    simulation.mesh.Lmat2 = alpha * identity;

    Field3d curl_b(simulation.domain.size(), 3);
    curl_b = simulation.mesh.curlB * simulation.fieldB;
    check::near(curl_b.norm(), 0.0, 1e-14, "constant B must have zero discrete curl");

    simulation.predict_electric_field(simulation.fieldEp, simulation.fieldE, simulation.fieldB, simulation.fieldJp);

    const auto sizes = simulation.fieldEp.sizes();
    for (int i = 0; i < sizes[0]; ++i)
        for (int j = 0; j < sizes[1]; ++j)
            for (int k = 0; k < sizes[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const double expected =
                        (e_value[d] - 0.5 * simulation.dt_ * j_value[d] + 0.5 * simulation.dt_ * curl_b(i, j, k, d)) /
                        (1.0 + alpha);
                    check::near(simulation.fieldEp(i, j, k, d), expected, 1e-11,
                                "diagonal Lmat2 must scale the predictor by 1/(1+alpha)");
                }

    const Field3d rhs = simulation.fieldE - 0.5 * simulation.dt_ * simulation.fieldJp + 0.5 * simulation.dt_ * curl_b;
    const Operator A = simulation.mesh.IMmat + simulation.mesh.Lmat2;
    const double residual = (A * simulation.fieldEp.data() - rhs.data()).norm() / rhs.data().norm();
    check::near(residual, 0.0, 1e-11, "diagonal-Lmat2 solve must satisfy (IMmat+Lmat2)*Ep = rhs");
}

void test_external_field_path() {
    char program[] = "field_solve";
    char* argv[] = {program};
    SimulationEcsim simulation(settings(), {{"particles", nlohmann::json::array()}}, 1, argv);
    simulation.init();

    const Vector3R e_value(0.7, -0.2, 0.4);
    const Vector3R b_value(0.3, -0.1, 0.5);
    const Vector3R j_value(0.4, 0.1, -0.3);
    const Vector3R e_external_value(0.2, -0.1, 0.05);
    fill_constant(simulation.fieldE, e_value);
    fill_constant(simulation.fieldB, b_value);
    fill_constant(simulation.fieldJp, j_value);
    fill_constant(simulation.fieldE_external, e_external_value);
    simulation.fieldEp.setZero();

    Operator identity(simulation.mesh.IMmat.rows(), simulation.mesh.IMmat.cols());
    identity.setIdentity();
    const double alpha = 0.25;
    simulation.mesh.Lmat2 = alpha * identity;

    Field3d curl_b(simulation.domain.size(), 3);
    curl_b = simulation.mesh.curlB * simulation.fieldB;
    check::near(curl_b.norm(), 0.0, 1e-14, "constant B must have zero discrete curl");

    // This is the external-field predictor used by make_step: Lmat2 carries
    // the correction and E_external enters the RHS through -Lmat2*E_external.
    simulation.predict_electric_field(simulation.fieldEp, simulation.fieldE, simulation.fieldE_external,
                                      simulation.fieldB, simulation.fieldJp);

    const auto sizes = simulation.fieldEp.sizes();
    for (int i = 0; i < sizes[0]; ++i)
        for (int j = 0; j < sizes[1]; ++j)
            for (int k = 0; k < sizes[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const double expected = (e_value[d] + 0.5 * simulation.dt_ * (curl_b(i, j, k, d) - j_value[d]) -
                                             alpha * e_external_value[d]) /
                                            (1.0 + alpha);
                    check::near(simulation.fieldEp(i, j, k, d), expected, 1e-11,
                                "external-field predictor must subtract Lmat2*E_external");
                }

    Field3d rhs = simulation.fieldE + 0.5 * simulation.dt_ * (curl_b - simulation.fieldJp) -
                  simulation.mesh.Lmat2 * simulation.fieldE_external;
    const Operator A = simulation.mesh.IMmat + simulation.mesh.Lmat2;
    const double residual = (A * simulation.fieldEp.data() - rhs.data()).norm() / rhs.data().norm();
    check::near(residual, 0.0, 1e-11, "external-field predictor must satisfy A*Ep = rhs");
}

void test_ecsim_corr_correct_e_constant_fields() {
    char program[] = "field_solve";
    char* argv[] = {program};
    SimulationEcsimCorr simulation(settings_corr(), {{"particles", nlohmann::json::array()}}, 1, argv);
    simulation.init();

    const Vector3R e_value(0.7, -0.2, 0.4);
    const Vector3R b_value(0.3, -0.1, 0.5);
    const Vector3R j_value(0.4, 0.1, -0.3);
    fill_constant(simulation.fieldE, e_value);
    fill_constant(simulation.fieldB, b_value);
    fill_constant(simulation.fieldJe, j_value);
    simulation.fieldEn.setZero();

    // The Corr solve_electric_field_from_deposited_current RHS is E + dt*(curlB*B - J) + Mmat*E. For constant
    // fields curlB*B=0 and Mmat*E=0, so the correction reduces to E - dt*J.
    const Field3d curl_b = simulation.mesh.curlB * simulation.fieldB;
    check::near(curl_b.norm(), 0.0, 1e-14, "constant B must have zero discrete curl");
    const Field3d m_times_e = simulation.mesh.Mmat * simulation.fieldE;
    check::near(m_times_e.norm(), 0.0, 1e-14, "constant E must be annihilated by Mmat");

    simulation.solve_electric_field_from_deposited_current(simulation.fieldEn, simulation.fieldE, simulation.fieldB,
                                                           simulation.fieldJe, simulation.dt_);

    const auto sizes = simulation.fieldEn.sizes();
    for (int i = 0; i < sizes[0]; ++i)
        for (int j = 0; j < sizes[1]; ++j)
            for (int k = 0; k < sizes[2]; ++k)
                for (int d = 0; d < 3; ++d) {
                    const double expected = e_value[d] - simulation.dt_ * j_value[d];
                    check::near(
                        simulation.fieldEn(i, j, k, d), expected, 1e-11,
                        "Corr solve_electric_field_from_deposited_current must give En = E - dt*J for constant fields");
                }

    const Field3d rhs = simulation.fieldE - simulation.dt_ * simulation.fieldJe;
    const double residual = (simulation.mesh.IMmat * simulation.fieldEn.data() - rhs.data()).norm() / rhs.data().norm();
    check::near(residual, 0.0, 1e-11,
                "Corr solve_electric_field_from_deposited_current must satisfy IMmat*En = E - dt*J");
}

}   // namespace

int main() {
    return check::run([] {
        check::scenario("constant fields with zero Lmat2", test_constant_fields_with_zero_lmat2);
        check::scenario("synthetic diagonal Lmat2 scales the solve", test_synthetic_diagonal_lmat2);
        check::scenario("external-field predictor subtracts Lmat2*E_external", test_external_field_path);
        check::scenario("Corr solve_electric_field_from_deposited_current gives En = E - dt*J for constant fields",
                        test_ecsim_corr_correct_e_constant_fields);
    });
}
