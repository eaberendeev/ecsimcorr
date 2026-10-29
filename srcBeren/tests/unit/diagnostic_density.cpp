#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>

#include "DensityMeasurements.h"
#include "DiagnosticMeasurements.h"
#include "OutputFactory.h"
#include "output_util.h"
#include "simulation_ecsim_corr.h"
#include "test_support.h"

namespace fs = std::filesystem;
namespace {
void requested_density_has_the_scheme_shape_and_current_particles() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(1, 1, 1));
    BoundaryConditionHandler boundaries;
    Species species;
    species.emplace("Beam", std::make_unique<ParticlesArray>(check::species("Beam", -2, 1, 1, 1), domain));
    species.emplace("Gas", std::make_unique<ParticlesArray>(check::species("Gas", 0, 1, 1, 1), domain));
    Field3d field(domain.size(), 3);
    field.setZero();
    for (auto shape : {ShapeType::Linear, ShapeType::Quadratic}) {
        DiagnosticMeasurements measured({domain, species, field, field, field, field, boundaries, shape});
        measured.request_density("Beam");
        measured.request_density("Beam");
        measured.request_density("Gas");
        // Requests do not sample early: newly created particles must be visible.
        for (const auto& [name, particles] : species)
            particles->add_particle(Particle(Vector3R(2.3, 2.4, 2.2), Vector3R(0, 0, 0)));
        measured.calculate();
        const int count = species.at("Beam")->get_total_num_of_particles();
        check::near(measured.density("Beam").data().sum(), -2 * count, 1e-13, "requested charge density");
        check::near(measured.density("Gas").data().sum(), count, 1e-13, "neutral number density");
        int support = 0;
        for (int i = 0; i < measured.density("Beam").capacity(); ++i) support += measured.density("Beam")[i] != 0;
        check::near(support, shape == ShapeType::Linear ? 8 : 27, 0, "scheme density shape support");
        const auto total = measure_charge_density(species, domain, boundaries, shape);
        check::near((total - measured.density("Beam")).norm(), 0, 1e-13, "neutral density entered electric charge");
    }
    bool unsupported = false;
    try {
        OutputFactory::create({{"outputs", nlohmann::json::array({{{"type", "charge_continuity"}, {"interval", 2}}})}},
                              domain, species, {{"Dt", 1}, {"Scheme", "ecsim"}}, field, field, field, field);
    } catch (const std::runtime_error& error) {
        unsupported = std::string(error.what()).find("trajectory current") != std::string::npos;
    }
    check::require(unsupported, "point-current ECSIM must reject the trajectory continuity diagnostic");
    DiagnosticMeasurements unused({domain, species, field, field, field, field, boundaries});
    unused.calculate();
    bool absent = false;
    try {
        unused.density("Beam");
    } catch (const std::out_of_range&) {
        absent = true;
    }
    check::require(absent, "unrequested density must not be deposited or cached");
}

nlohmann::json settings() {
    return {{"Scheme", "ecsim_corr"},
            {"Dt", 0.25},
            {"NumCellsX", 6},
            {"NumCellsY", 6},
            {"NumCellsZ", 6},
            {"Dx", 1},
            {"Dy", 1},
            {"Dz", 1},
            {"StartFromTime", 0},
            {"StartTimeStep", 0},
            {"LastTimestep", 5},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collider", "None"},
            {"Collisions", nlohmann::json::array()},
            {"Boundary_conditions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::parse(R"([
                {"type":"charge_continuity","interval":2},
                {"type":"full3d_all","timesteps":[2,4]},
                {"type":"recovery","interval":4}
            ])")}}}};
}

void check_continuity_file(const std::vector<double>& times) {
    std::ifstream file("continuity.txt");
    std::string header;
    std::getline(file, header);
    check::require(header == "Time L2 Linf", "continuity column format");
    for (double expected_time : times) {
        double time, l2, linf;
        check::require(bool(file >> time >> l2 >> linf), "missing selected continuity row");
        check::near(time, expected_time, 0, "continuity cadence");
        check::near(l2, 0, 1e-13, "stationary-particle continuity L2");
        check::near(linf, 0, 1e-13, "stationary-particle continuity Linf");
    }
    file >> std::ws;
    check::require(file.peek() == std::char_traits<char>::eof(), "continuity wrote step 0 or an unchecked step");
}

void selected_steps_capture_both_endpoints_and_restart(const fs::path& root) {
    fs::current_path(root);
    auto config = settings();
    const auto particles = nlohmann::json{{"particles", nlohmann::json::array({check::species("Beam", -1, 1, 1, 1)})}};
    char program[] = "diagnostic_density";
    char* argv[] = {program};
    auto sim = build_simulation(config, particles, 1, argv);
    sim->init();
    sim->make_diagnostic(0);
    auto& beam = *sim->species.at("Beam");
    auto& corr = dynamic_cast<SimulationEcsimCorr&>(*sim);
    for (int step = 1; step <= 5; ++step) {
        sim->prepare_step(step);
        // A source before transport belongs in BOTH endpoints of this step.
        beam.add_particle(Particle(Vector3R(2.3, 2.4, 2.2), Vector3R(0, 0, 0)));
        sim->make_step(step);
        if (step % 2 == 0)
            // A particle created after transport belongs to output, not to the
            // completed transport residual (as with emission/product creation).
            beam.add_particle(Particle(Vector3R(2.3, 2.4, 2.2), Vector3R(0, 0, 0)));
        Diagnostics values(sim->species);
        corr.diagnostic_energy(values);
        check::require(values.charge_continuity.has_value() == (step % 2 == 0),
                       "continuity result survived a skipped step or missed a selected step");
        if (values.charge_continuity)
            check::near(values.charge_continuity->l2, 0, 1e-13, "endpoints span a diagnostic gap or include a source");
        sim->make_diagnostic(step);
        if (step % 2 == 0) {
            Field3d density(sim->domain.size(), 1);
            read_field_from_file("Full3D/Particles/Beam/Density_" + std::to_string(step) + ".3d", density);
            check::near(density.data().sum(), -beam.get_total_num_of_particles(), 1e-6,
                        "output density must include particles created after transport");
        }
    }
    sim->finalize();
    sim->finalize();
    check_continuity_file({0.5, 1});
    sim.reset();

    // Recovery contains particles and fields, not a previous density snapshot.
    fs::create_directory(root / "resumed");
    fs::current_path(root / "resumed");
    config["StartFromTime"] = 4;
    config["StartTimeStep"] = 4;
    config["LastTimestep"] = 6;
    auto resumed = build_simulation(config, particles, 1, argv);
    resumed->init();
    resumed->calculate();
    resumed->finalize();
    check_continuity_file({1.5});
}

void continuity_rejects_nonfinite_samples() {
    Field3d before(Vector3I(2, 2, 2), 1), after(Vector3I(2, 2, 2), 1), current(Vector3I(2, 2, 2), 3);
    before.setZero();
    after.data().setConstant(2);
    current.setZero();
    Operator divergence(before.capacity(), current.capacity());
    const auto result = measure_charge_continuity(before, after, current, divergence, 0.5);
    check::near(result.linf, 4, 0, "continuity uses the full single-step dt");
    check::near(result.l2, 4 * std::sqrt(8), 1e-13, "continuity full-grid L2 norm");
    after[0] = std::numeric_limits<double>::quiet_NaN();
    bool rejected = false;
    try {
        measure_charge_continuity(before, after, current, divergence, 0.5);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    check::require(rejected, "nonfinite continuity must fail");
}

void rectangular_divergence_and_nonfinite_reduction() {
    Field3d before(Vector3I(2, 2, 2), 1), after(Vector3I(2, 2, 2), 1), current(Vector3I(2, 2, 2), 3);
    before.setZero();
    current.setZero();
    Operator divergence(before.capacity(), current.capacity());
    for (int row = 0; row < before.capacity(); ++row) {
        divergence.insert(row, 3 * row) = 2;
        divergence.insert(row, 3 * row + 1) = -1;
        current[3 * row] = row;
        current[3 * row + 1] = row + 1;
    }
    divergence.makeCompressed();
    const int original_threads = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        after.data().setConstant(2);
        // drho/dt=4 and div J=row-1, hence the eight residuals are 3,...,10.
        const ChargeContinuityResidual result = measure_charge_continuity(before, after, current, divergence, 0.5);
        check::near(result.linf, 10, 0, "rectangular divergence must produce a scalar residual at each node");
        check::near(result.l2, std::sqrt(380), 1e-13, "independent norm of the eight continuity residuals");
        for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()}) {
            // The last node belongs to another worker when multiple threads run.
            after[after.capacity() - 1] = bad;
            bool rejected = false;
            try {
                measure_charge_continuity(before, after, current, divergence, 0.5);
            } catch (const std::runtime_error&) {
                rejected = true;
            }
            check::require(rejected, "parallel continuity scan missed a nonfinite value on a worker");
        }
    }
    omp_set_num_threads(original_threads);
}

}   // namespace

int main() {
    const auto original = fs::current_path();
    const auto root =
        fs::temp_directory_path() /
        ("beren-density-test-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directory(root);
    const int result = check::run([&] {
        check::scenario("Requested densities: shape, charge/number and current particle state",
                        requested_density_has_the_scheme_shape_and_current_particles);
        check::scenario("Selected continuity steps: endpoint timing, output density and recovery",
                        [&] { selected_steps_capture_both_endpoints_and_restart(root); });
        check::scenario("Continuity norms use one-step dt and reject NaN", continuity_rejects_nonfinite_samples);
        check::scenario("Rectangular divergence and NaN/Inf rejection with 1/2/4 threads",
                        rectangular_divergence_and_nonfinite_reduction);
    });
    fs::current_path(original);
    fs::remove_all(root);
    return result;
}
