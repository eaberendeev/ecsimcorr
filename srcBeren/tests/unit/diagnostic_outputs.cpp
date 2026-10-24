#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "DiagnosticOutput.h"

namespace fs = std::filesystem;
namespace {
int checks = 0;
void require(bool condition, const std::string& message) {
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
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
    require(width > 0 && height > 0 && !result.empty(), "invalid plane header/data");
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

void exercise(const std::string& type, int mode, const fs::path& root) {
    const bool ordinary = (mode & 1) != 0;
    const bool averaged = (mode & 2) != 0;
    const fs::path work = root / (type + std::to_string(mode));
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
    Diagnostics diagnostic(nlohmann::json::object(), domain, species);
    Field3d field(domain.size(), 3), zero(domain.size(), 3);
    zero.setZero();
    auto output_cfg = nlohmann::json::parse(R"({"fields":["En"],"species":["Electrons"],"planes":{"z":[3]}})");
    output_cfg["type"] = type;
    if (ordinary)
        output_cfg["interval"] = 2;
    if (averaged)
        output_cfg["average_steps"] = 2;
    if (mode == 4)
        output_cfg["average_steps"] = 0;
    auto outputs = OutputFactory::create({{"outputs", nlohmann::json::array({output_cfg})}}, domain, species,
                                         {{"Dt", 1.0}}, field, zero, zero, zero);
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
            type + ": unexpected files in mode " + std::to_string(mode));
    for (const auto& file : files) {
        const auto name = file.filename().string();
        require((name.find("_avg") != std::string::npos) ? averaged : ordinary, "unrequested output");
        const auto split = name.find("_Plane");
        const std::string basename = name.substr(0, split);
        const fs::path reference_dir = root / (type + "1") / fs::relative(file.parent_path(), work);
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
}   // namespace

int main() {
    const auto original = fs::current_path();
    const auto root =
        fs::temp_directory_path() /
        ("beren-output-test-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directory(root);
    int result = 0;
    try {
        for (const std::string type : {"fields_2d", "density_2d", "current_2d", "pressure_2d"}) {
            // Reference instantaneous files must exist before average-only checks.
            for (int mode : {1, 0, 2, 3, 4}) exercise(type, mode, root);
            Domain domain;
            domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
            Species species;
            Field3d field(domain.size(), 3);
            for (const auto& invalid : {nlohmann::json{{"interval", 0}}, nlohmann::json{{"interval", -1}},
                                        nlohmann::json{{"average_steps", -1}}}) {
                auto entry = invalid;
                entry["type"] = type;
                bool rejected = false;
                try {
                    OutputFactory::create({{"outputs", nlohmann::json::array({entry})}}, domain, species, {{"Dt", 1}},
                                          field, field, field, field);
                } catch (const std::runtime_error&) {
                    rejected = true;
                }
                require(rejected, type + ": invalid settings accepted");
            }
        }
        Domain domain;
        Species species;
        Field3d field;
        require(
            OutputFactory::create(nlohmann::json::object(), domain, species, {{"Dt", 1}}, field, field, field, field)
                .empty(),
            "absent outputs should disable diagnostics");
        require(OutputFactory::create({{"outputs", nlohmann::json::array()}}, domain, species, {{"Dt", 1}}, field,
                                      field, field, field)
                    .empty(),
                "empty outputs should disable diagnostics");
        std::cout << checks << " output checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        result = 1;
    }
    fs::current_path(original);
    fs::remove_all(root);
    return result;
}
