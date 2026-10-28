#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "Diagnostic.h"
#include "OutputFactory.h"
#include "diagnostic_output_validation.h"
#include "output_util.h"
#include "test_case.h"

namespace fs = std::filesystem;
namespace {
using check::require;

enum class OutputMode { Snapshots, Disabled, Averages, Both, ZeroAverageWindow };
const char* mode_name(OutputMode mode) {
    switch (mode) {
        case OutputMode::Snapshots:
            return "snapshots";
        case OutputMode::Disabled:
            return "disabled";
        case OutputMode::Averages:
            return "averages";
        case OutputMode::Both:
            return "snapshots_and_averages";
        case OutputMode::ZeroAverageWindow:
            return "zero_average_window";
    }
    throw std::runtime_error("unknown output mode");
}

std::vector<float> read_plane(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    require(input.is_open(), "missing plane: " + path.string());
    float width = 0, height = 0;
    input.read(reinterpret_cast<char*>(&width), sizeof(float));
    input.read(reinterpret_cast<char*>(&height), sizeof(float));
    std::vector<float> result;
    float value;
    while (input.read(reinterpret_cast<char*>(&value), sizeof(float))) result.push_back(value);
    const int side = 4 + 2 * GHOST_CELLS + 1;
    require(width == side && height == side, "wrong plane dimensions");
    require(input.gcount() == 0 && result.size() == size_t(3 * side * side), "truncated or oversized plane payload");
    return result;
}

void check_scaled(const fs::path& path, const std::vector<float>& reference, double factor) {
    const auto actual = read_plane(path);
    require(actual.size() == reference.size(), "plane shape changed");
    for (size_t i = 0; i < actual.size(); ++i)
        require(std::isfinite(actual[i]) &&
                    std::abs(actual[i] - factor * reference[i]) < 1e-5 * (1 + std::abs(factor * reference[i])),
                "wrong mean: " + path.string());
}

void verify_output_files_and_window_means(const std::string& type, OutputMode mode, const fs::path& root) {
    const bool ordinary = mode == OutputMode::Snapshots || mode == OutputMode::Both;
    const bool averaged = mode == OutputMode::Averages || mode == OutputMode::Both;
    const fs::path work = root / (type + "_" + mode_name(mode));
    fs::create_directory(work);
    fs::current_path(work);
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    Species species;
    auto cfg = nlohmann::json::parse(
        R"({"Name":"Electrons","Charge":-1,"Mass":1,"Density":1,"NumPartPerCell":1,"distribution":[]})");
    species.emplace("Electrons", std::make_unique<ParticlesArray>(cfg, domain));
    auto& sp = *species.at("Electrons");
    sp.add_particle(Particle(Vector3R(2.3, 2.4, 2.0), Vector3R(1, 2, 3)));
    Diagnostics diagnostic(species);
    Field3d field(domain.size(), 3), zero(domain.size(), 3);
    zero.setZero();
    auto output_cfg = nlohmann::json::parse(R"({"fields":["En"],"species":["Electrons"],"planes":{"z":[3]}})");
    output_cfg["type"] = type;
    if (ordinary)
        output_cfg["interval"] = 2;
    if (averaged)
        output_cfg["average_steps"] = 2;
    if (mode == OutputMode::ZeroAverageWindow)
        output_cfg["average_steps"] = 0;
    auto outputs = OutputFactory::create({{"outputs", nlohmann::json::array({output_cfg})}}, domain, species,
                                         {{"Dt", 1.0}}, field, zero, zero, zero);
    // Samples 1,3,5: the full two-sample mean is 2, the residual mean is 5.
    // Scaling v by sqrt(sample) scales the raw second moment by sample too.
    for (int step = 0; step < 3; ++step) {
        const double value = 2 * step + 1;
        field.data().setConstant(value);
        sp.densityOnGrid.data().setConstant(value);
        sp.currentOnGrid.data().setConstant(value);
        for (int cell = 0; cell < sp.size(); ++cell)
            for (auto& particle : sp.particlesData(cell)) particle.velocity = std::sqrt(value) * Vector3R(1, 2, 3);
        const Field3d field_before = field, density_before = sp.densityOnGrid, current_before = sp.currentOnGrid;
        const auto values_before = diagnostic.energy;
        const double kinetic_before = sp.get_kinetic_energy();
        const int count_before = sp.get_total_num_of_particles();
        const Diagnostics& measured = diagnostic;
        outputs.front()->output(step, measured);
        check::near((field - field_before).norm(), 0, 0, "output changed the electric field");
        check::near((sp.densityOnGrid - density_before).norm(), 0, 0, "output changed the deposited density");
        check::near((sp.currentOnGrid - current_before).norm(), 0, 0, "output changed the deposited current");
        check::near(sp.get_kinetic_energy(), kinetic_before, 0, "output changed particle kinetic energy");
        require(sp.get_total_num_of_particles() == count_before, "output changed the particle count");
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const auto& particle : sp.particlesData(cell))
                for (int axis = 0; axis < 3; ++axis) {
                    check::near(particle.coord[axis], Vector3R(2.3, 2.4, 2)[axis], 0, "output moved a particle");
                    check::near(particle.velocity[axis], std::sqrt(value) * Vector3R(1, 2, 3)[axis], 0,
                                "output changed a particle velocity");
                }
        require(diagnostic.energy == values_before, "output changed measured diagnostic values");
    }
    outputs.front()->finalize();
    outputs.front()->finalize();
    std::vector<fs::path> files;
    for (const auto& entry : fs::recursive_directory_iterator(work))
        if (entry.is_regular_file())
            files.push_back(entry.path());
    const int components = type == "pressure_2d" ? 6 : 1;
    require(files.size() == static_cast<size_t>(components * (2 * ordinary + 2 * averaged)),
            type + ": unexpected files in mode " + mode_name(mode));
    for (const auto& file : files) {
        const auto name = file.filename().string();
        require((name.find("_avg") != std::string::npos) ? averaged : ordinary, "unrequested output");
        const auto split = name.find("_Plane");
        const std::string basename = name.substr(0, split);
        const fs::path reference_dir =
            root / (type + "_" + mode_name(OutputMode::Snapshots)) / fs::relative(file.parent_path(), work);
        const auto reference = read_plane(reference_dir / (basename + "_PlaneZ_003_0000"));
        const double factor = name.ends_with("avg0000") ? 2.0 : name.ends_with("0001") ? 5.0 : 1.0;
        check_scaled(file, reference, factor);
        if (type == "pressure_2d" && basename == "Pzz") {
            bool nonzero = false;
            for (float value : reference) nonzero |= value > 0;
            require(nonzero, "pressure test must contain nonzero moment");
        }
    }
}

void absent_output_list_disables_outputs(const fs::path& root) {
    const auto work = root / "outputs_disabled";
    fs::create_directory(work);
    fs::current_path(work);
    Domain domain;
    Species species;
    Field3d field;
    Diagnostics diagnostic(species);
    require(fs::is_empty(work), "constructing diagnostic values must not create output directories");
    require(OutputFactory::create(nlohmann::json::object(), domain, species, {{"Dt", 1}}, field, field, field, field)
                .empty(),
            "absent outputs should disable diagnostics");
    require(OutputFactory::create({{"outputs", nlohmann::json::array()}}, domain, species, {{"Dt", 1}}, field, field,
                                  field, field)
                .empty(),
            "empty outputs should disable diagnostics");
    require(fs::is_empty(work), "disabled outputs must not create files or directories");
}

void tables_and_console_preserve_prepared_values(const fs::path& root) {
    const auto work = root / "tables_and_console";
    fs::create_directory(work);
    fs::current_path(work);
    Domain domain;
    Species species;
    Diagnostics diagnostic(species);
    diagnostic.addEnergy("known", 2.5);
    diagnostic.addBoundary("lost", 0.75);
    // Missing values retain the legacy zero in the files; reading must not
    // insert them or the console's optional energy fields into the measurement.
    diagnostic.energyOrder.push_back("missing");
    diagnostic.boundaryOrder.push_back("missing_boundary");
    const auto values = diagnostic.energy;
    const auto energy_columns = diagnostic.energyOrder, boundary_columns = diagnostic.boundaryOrder;
    Field3d field;
    const auto settings = nlohmann::json::parse(
        R"({"outputs":[{"type":"energy_balance"},{"type":"boundary_stats"},{"type":"console_summary"}]})");
    auto outputs = OutputFactory::create(settings, domain, species, {{"Dt", 0.5}}, field, field, field, field);
    const Diagnostics& measured = diagnostic;
    for (int step : {0, 1})
        for (auto& output : outputs) output->output(step, measured);
    for (auto& output : outputs) {
        output->finalize();
        output->finalize();
    }
    require(diagnostic.energy == values, "table/console output changed the diagnostic values");
    require(diagnostic.energyOrder == energy_columns && diagnostic.boundaryOrder == boundary_columns,
            "table/console output changed the column schema");
    for (const std::string filename : {"energy.txt", "boundary.txt"}) {
        std::ifstream input(filename);
        std::string header;
        std::getline(input, header);
        require(header == (filename == "energy.txt" ? "Time known missing " : "Time lost missing_boundary "),
                "table header or column order changed");
        for (int step : {0, 1}) {
            double time = -1, known = -1, missing = -1;
            require(bool(input >> time >> known >> missing), "table row missing");
            check::near(time, step * 0.5, 0, "table time normalization");
            check::near(known, filename == "energy.txt" ? 2.5 : 0.75, 0, "prepared table value");
            check::near(missing, 0, 0, "legacy missing table value");
        }
        input >> std::ws;
        require(input.peek() == std::char_traits<char>::eof(), "extra table rows or values");
    }
}

void configured_dumps_preserve_values_and_cadence(const fs::path& root) {
    const auto work = root / "scheduled_dumps";
    fs::create_directory(work);
    fs::current_path(work);
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    Species species;
    const auto config = nlohmann::json::parse(
        R"({"Name":"Electrons","Charge":-1,"Mass":1,"Density":1,"NumPartPerCell":1,"distribution":[]})");
    species.emplace("Electrons", std::make_unique<ParticlesArray>(config, domain));
    auto& electrons = *species.at("Electrons");
    electrons.add_particle(Particle(Vector3R(2.3, 2.4, 2), Vector3R(1, 2, 3)));
    electrons.densityOnGrid.data().setConstant(-2);
    electrons.currentOnGrid.data().setConstant(0.75);
    Diagnostics diagnostic(species);
    Field3d en(domain.size(), 3), bn(domain.size(), 3), external_e(domain.size(), 3), external_b(domain.size(), 3);
    en.data().setConstant(2);
    bn.data().setConstant(3);
    external_e.data().setConstant(4);
    external_b.data().setConstant(5);
    const auto settings = nlohmann::json::parse(R"({"outputs":[
        {"type":"fields_2d","planes":{"z":[3]},"interval":2},
        {"type":"fields_3d","fields":["E","Bn","E_ex"],"timesteps":[2]},
        {"type":"full3d_all","timesteps":[2]},
        {"type":"charge_density_2d","planes":{"z":[3]},"interval":2},
        {"type":"total_current_2d","planes":{"z":[3]},"interval":2},
        {"type":"energy_spectrum","species":["Electrons"],"interval":2},
        {"type":"recovery","interval":2},
        {"type":"probes","points":[[2,2,2]],"species":["Electrons"],"interval":2}
    ]})");
    auto outputs = OutputFactory::create(settings, domain, species, {{"Dt", 1}}, en, bn, external_e, external_b);
    for (auto& output : outputs) output->output(1, diagnostic);
    // Probe preparation opens an empty file; no output may write a sample at step 1.
    for (const auto& entry : fs::recursive_directory_iterator(work))
        if (entry.is_regular_file())
            require(fs::file_size(entry.path()) == 0, "dump ignored its configured cadence");
    for (auto& output : outputs) output->output(2, diagnostic);
    for (auto& output : outputs) {
        output->finalize();
        output->finalize();
    }

    const auto verify_dump = [&](const std::string& path, const Field3d& expected) {
        Field3d actual(expected.sizes(), expected.nd());
        read_field_from_file(path, actual);
        check::near((actual - expected).norm(), 0, 0, "wrong field in " + path);
    };
    const Field3d total_e = en + external_e, total_b = bn + external_b;
    verify_dump("Fields/Diag3D/FieldE_2.3d", total_e);
    verify_dump("Fields/Diag3D/FieldBn_2.3d", bn);
    verify_dump("Fields/Diag3D/FieldE_external_2.3d", external_e);
    verify_dump("Full3D/Fields/E_2.3d", total_e);
    verify_dump("Full3D/Fields/B_2.3d", total_b);
    verify_dump("Full3D/Fields/En_2.3d", en);
    verify_dump("Full3D/Fields/Bn_2.3d", bn);
    verify_dump("Full3D/Fields/E_ex_2.3d", external_e);
    verify_dump("Full3D/Fields/B_init_2.3d", external_b);
    verify_dump("Full3D/Particles/Electrons/Density_2.3d", electrons.densityOnGrid);
    verify_dump("Full3D/Particles/Electrons/Current_2.3d", electrons.currentOnGrid);
    verify_dump("Recovery/Fields/FieldE.backup", en);
    verify_dump("Recovery/Fields/FieldB.backup", bn);
    std::ifstream particles("Recovery/Particles/Electrons/Electrons.backup", std::ios::binary);
    int count = 0;
    Particle stored;
    particles.read(reinterpret_cast<char*>(&count), sizeof(count));
    particles.read(reinterpret_cast<char*>(&stored), sizeof(stored));
    require(bool(particles) && count == 1, "recovery must contain the initial particle");
    for (int axis = 0; axis < 3; ++axis) {
        check::near(stored.coord[axis], Vector3R(2.3, 2.4, 2)[axis], 0, "recovery particle position");
        check::near(stored.velocity[axis], Vector3R(1, 2, 3)[axis], 0, "recovery particle velocity");
    }

    for (float value : read_plane("Fields/Diag2D/FieldE_PlaneZ_003_0001"))
        check::near(value, 6, 0, "default E slice includes the external field");
    for (float value : read_plane("Fields/Diag2D/FieldB_PlaneZ_003_0001"))
        check::near(value, 8, 0, "default B slice includes the prescribed field");
    const auto charge = read_plane("Fields/Diag2D/ChargeDensity_PlaneZ_003_0001");
    const auto current = read_plane("Fields/Diag2D/TotalCurrent_PlaneZ_003_0001");
    const size_t area = charge.size() / 3;
    for (size_t i = 0; i < charge.size(); ++i) {
        check::near(charge[i], i < area ? -2 : 0, 0, "total charge slice and zero unused components");
        check::near(current[i], 0.75, 0, "total current slice");
    }
    std::ifstream probes("Probes/probes.txt");
    std::string header;
    std::getline(probes, header);
    double value = 0;
    const std::vector<double> expected_probe{2, 2, 2, 2, 6, 6, 6, 8, 8, 8, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, -2};
    for (double expected : expected_probe) {
        require(bool(probes >> value), "missing probe value");
        check::near(value, expected, 1e-12, "uniform fields and density at the probe point");
    }
    probes >> std::ws;
    require(probes.peek() == std::char_traits<char>::eof(), "unexpected probe sample");

    std::ifstream spectrum("Particles/Electrons/energy_spectrum_000002.txt");
    double energy = 0, bin_count = 0, total_count = 0;
    int bins = 0;
    while (spectrum >> energy >> bin_count) {
        require(std::isfinite(energy) && std::isfinite(bin_count), "nonfinite energy spectrum");
        if (bin_count > 0)
            check::near(energy, 7, 2e-4, "particle energy in occupied spectrum bin");
        total_count += bin_count;
        ++bins;
    }
    check::near(total_count, 1, 0, "spectrum must count the initial particle once");
    require(bins == 1000, "spectrum bin count changed");
    check::near(electrons.get_kinetic_energy(), 7, 0, "scheduled dumps changed particle energy");
}

}   // namespace

int main() {
    const auto original = fs::current_path();
    const auto root =
        fs::temp_directory_path() /
        ("beren-output-test-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directory(root);
    const int result = check::run([&] {
        for (const std::string type : {"fields_2d", "density_2d", "current_2d", "pressure_2d"}) {
            // Instantaneous reference files must exist before averaged checks.
            for (auto mode : {OutputMode::Snapshots, OutputMode::Disabled, OutputMode::Averages, OutputMode::Both,
                              OutputMode::ZeroAverageWindow}) {
                check::scenario(type + ": " + mode_name(mode) + ", step 0 and full/residual windows",
                                [&] { verify_output_files_and_window_means(type, mode, root); });
            }
        }
        verify_diagnostic_output_failures(root);
        check::scenario("Absent/empty outputs and diagnostic values create no files or directories",
                        [&] { absent_output_list_disables_outputs(root); });
        check::scenario("Tables and console read prepared values without inserting missing entries",
                        [&] { tables_and_console_preserve_prepared_values(root); });
        check::scenario("Scheduled fields, charge/current, recovery, probes and spectrum preserve values and cadence",
                        [&] { configured_dumps_preserve_values_and_cadence(root); });
    });
    fs::current_path(original);
    fs::remove_all(root);
    return result;
}
