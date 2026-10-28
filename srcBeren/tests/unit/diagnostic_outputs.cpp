#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "DiagnosticOutput.h"
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
        outputs.front()->output(step, diagnostic);
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

void invalid_output_settings_are_rejected(const std::string& type) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    Species species;
    Field3d field(domain.size(), 3);
    for (const auto& invalid :
         {nlohmann::json{{"interval", 0}}, nlohmann::json{{"interval", -1}}, nlohmann::json{{"average_steps", -1}}}) {
        auto entry = invalid;
        entry["type"] = type;
        bool rejected = false;
        try {
            OutputFactory::create({{"outputs", nlohmann::json::array({entry})}}, domain, species, {{"Dt", 1}}, field,
                                  field, field, field);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        require(rejected, type + ": invalid settings accepted");
    }
}

void absent_output_list_disables_outputs() {
    Domain domain;
    Species species;
    Field3d field;
    require(OutputFactory::create(nlohmann::json::object(), domain, species, {{"Dt", 1}}, field, field, field, field)
                .empty(),
            "absent outputs should disable diagnostics");
    require(OutputFactory::create({{"outputs", nlohmann::json::array()}}, domain, species, {{"Dt", 1}}, field, field,
                                  field, field)
                .empty(),
            "empty outputs should disable diagnostics");
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
            check::scenario(type + ": invalid intervals and windows must fail",
                            [&] { invalid_output_settings_are_rejected(type); });
        }
        check::scenario("Absent/empty output lists produce no diagnostics", absent_output_list_disables_outputs);
    });
    fs::current_path(original);
    fs::remove_all(root);
    return result;
}
