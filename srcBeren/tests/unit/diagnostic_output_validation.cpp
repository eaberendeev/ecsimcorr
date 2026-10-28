#include "diagnostic_output_validation.h"

#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "Diagnostic.h"
#include "DiagnosticMeasurements.h"
#include "OutputFactory.h"
#include "output_util.h"
#include "outputs/SliceFiles.h"
#include "recovery.h"
#include "test_support.h"

namespace fs = std::filesystem;
namespace {
void failed_writes_are_reported(const fs::path& root, const Domain& domain, const Species& species,
                                const Field3d& field) {
    const auto work = root / "failed_writes";
    fs::create_directory(work);
    fs::current_path(work);
    Diagnostics measured(species);
    const auto expect_failure = [](const auto& write, const std::string& filename) {
        bool failed = false;
        try {
            write();
        } catch (const std::runtime_error& error) {
            failed = std::string(error.what()).find(filename) != std::string::npos;
        }
        check::require(failed, "write failure must identify " + filename);
    };
    // /dev/full opens successfully but rejects writes, including buffered flushes.
    for (const auto& [type, filename] : std::vector<std::pair<std::string, fs::path>>{
             {"energy_balance", "energy.txt"},
             {"boundary_stats", "boundary.txt"},
             {"probes", "Probes/probes.txt"},
             {"energy_spectrum", "Particles/Electrons/energy_spectrum_000000.txt"},
             {"recovery", "Recovery/Particles/Electrons/Electrons.backup"}}) {
        if (!filename.parent_path().empty())
            fs::create_directories(filename.parent_path());
        fs::create_symlink("/dev/full", filename);
        const auto settings = nlohmann::json::parse(R"({"interval":1,"points":[[1,1,1]]})");
        auto entry = settings;
        entry["type"] = type;
        auto outputs = OutputFactory::create({{"outputs", nlohmann::json::array({entry})}}, domain, species,
                                             {{"Dt", 1}}, field, field, field, field);
        BoundaryConditionHandler boundaries;
        expect_failure(
            [&] { write_diagnostics(0, measured, {domain, species, field, field, field, field, boundaries}, outputs); },
            filename.filename().string());
    }
    expect_failure([&] { output_array2d({1, 2, 3}, 1, 1, "/dev/full"); }, "/dev/full");
    expect_failure([&] { write_field_to_file("/dev/full", field); }, "/dev/full");
}

void plane_axes_preserve_coordinates(const fs::path& root) {
    const auto work = root / "plane_axes";
    fs::create_directory(work);
    fs::current_path(work);
    Field3d field(Vector3I(6, 7, 8), 3);
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 7; ++j)
            for (int k = 0; k < 8; ++k)
                for (int d = 0; d < 3; ++d) field(i, j, k, d) = 1000 * d + 100 * i + 10 * j + k;
    const Vector3I start(1, 2, 3), end(5, 5, 5);
    for (int axis = 0; axis < 3; ++axis) {
        const int size1 = axis == 0 ? 3 : 4;
        const int size2 = axis == 2 ? 3 : 2;
        output_field_plane(field, start, end, 1, axis, 3, "field", "0");
        std::ifstream file(std::string("field_Plane") + "XYZ"[axis] + "_001_0", std::ios::binary);
        float width = 0, height = 0;
        file.read(reinterpret_cast<char*>(&width), sizeof(float));
        file.read(reinterpret_cast<char*>(&height), sizeof(float));
        check::require(width == size1 && height == size2, "plane dimensions");
        for (int d = 0; d < 3; ++d)
            for (int u = 0; u < size1; ++u)
                for (int v = 0; v < size2; ++v) {
                    Vector3I node = start;
                    node[axis] += 1;
                    int coordinate = 0;
                    for (int a = 0; a < 3; ++a)
                        if (a != axis)
                            node[a] += coordinate++ == 0 ? u : v;
                    float actual = 0;
                    file.read(reinterpret_cast<char*>(&actual), sizeof(float));
                    check::require(bool(file), "truncated plane");
                    check::near(actual, 1000 * d + 100 * node.x() + 10 * node.y() + node.z(), 0,
                                "plane coordinate/component ordering");
                }
    }
}

struct InvalidSettings {
    std::string name;
    nlohmann::json config;
    std::string option;
};

void require_config_failure(const nlohmann::json& config, const std::string& option, const Domain& domain,
                            const Species& species, const Field3d& field) {
    bool rejected = false;
    try {
        OutputFactory::create(config, domain, species, {{"Dt", 1}, {"Scheme", "ecsim_corr"}}, field, field, field,
                              field);
    } catch (const std::runtime_error& error) {
        rejected = std::string(error.what()).find(option) != std::string::npos;
    }
    check::require(rejected, "invalid settings must fail at creation and identify '" + option + "'");
}

void invalid_settings_create_no_output(const fs::path& root, const Domain& domain, const Species& species,
                                       const Field3d& field) {
    using Json = nlohmann::json;
    std::vector<InvalidSettings> cases;
    const auto add = [&](const std::string& name, const Json& entry, const std::string& option) {
        cases.push_back({name, {{"outputs", Json::array({entry})}}, option});
    };
    for (const std::string type : {"fields_2d", "density_2d", "current_2d", "pressure_2d"}) {
        for (const auto& interval : {Json(0), Json(-1), Json(1.5), Json("2"), Json(true), Json(2147483648LL)})
            add(type + ": invalid snapshot interval " + interval.dump(), {{"type", type}, {"interval", interval}},
                "interval");
        for (const auto& window : {Json(-1), Json(1.5), Json("2"), Json(2147483648LL)})
            add(type + ": invalid averaging window " + window.dump(), {{"type", type}, {"average_steps", window}},
                "average_steps");
    }
    for (const auto& interval : {Json(0), Json(-1), Json(1.5), Json("2"), Json(true), Json(2147483648LL)})
        add("invalid continuity interval " + interval.dump(), {{"type", "charge_continuity"}, {"interval", interval}},
            "interval");
    for (const std::string type : {"density_2d", "current_2d", "pressure_2d", "energy_spectrum", "probes"}) {
        for (const auto& selection : {Json("Electrons"), Json(3), Json::array({"missing"}), Json::array({3})})
            add(type + ": invalid species selection " + selection.dump(), {{"type", type}, {"species", selection}},
                "species");
    }
    for (const std::string type : {"fields_2d", "fields_3d"}) {
        for (const auto& fields : {Json("E"), Json::array({"typo"}), Json::array({4})})
            add(type + ": invalid field selection " + fields.dump(), {{"type", type}, {"fields", fields}}, "fields");
    }
    add("3D output rejects a field it does not implement", {{"type", "fields_3d"}, {"fields", {"B_init"}}}, "fields");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const auto& planes : {Json::array(), Json{{"q", {1}}}, Json{{"z", 3}}, Json{{"z", {"3"}}}, Json{{"z", {nan}}},
                               Json{{"z", {inf}}}, Json{{"z", {1e300}}}})
        add("invalid plane coordinates " + planes.dump(), {{"type", "fields_2d"}, {"planes", planes}}, "planes");
    for (const auto& points : {Json(3), Json::array({{1, 2}}), Json::array({{1, 2, 3, 4}}), Json::array({{1, "x", 2}}),
                               Json::array({{nan, 2, 2}}), Json::array({{inf, 2, 2}}), Json::array({{-2, 2, 2}}),
                               Json::array({{100, 2, 2}})})
        add("invalid probe point " + points.dump(), {{"type", "probes"}, {"points", points}}, "points");
    for (const auto& steps : {Json(2), Json::array({-1}), Json::array({1.5}), Json::array({2147483648LL})})
        add("invalid output steps " + steps.dump(), {{"type", "fields_3d"}, {"timesteps", steps}}, "timesteps");
    for (const auto& entry : {Json::object(), Json(1), Json{{"type", 3}}, Json{{"type", "typo"}}})
        add("invalid output type " + entry.dump(), entry, "type");
    cases.push_back({"diagnostics must be an object", Json::array(), "diagnostics"});
    cases.push_back({"outputs must be an array", {{"outputs", 1}}, "outputs"});

    int index = 0;
    for (const auto& invalid : cases) {
        check::scenario(invalid.name, [&] {
            const auto work = root / ("invalid_settings_" + std::to_string(index++));
            fs::create_directory(work);
            fs::current_path(work);
            require_config_failure(invalid.config, invalid.option, domain, species, field);
            check::require(fs::is_empty(work), "invalid configuration opened output files or created directories");
        });
    }
}

void blocked_directories_fail_before_writing(const fs::path& root, const Domain& domain, const Species& species,
                                             const Field3d& field) {
    const std::vector<std::pair<std::string, fs::path>> cases{{"recovery", "Recovery/Particles/Electrons"},
                                                              {"full3d_all", "Full3D/Particles/Electrons"},
                                                              {"fields_3d", "Fields/Diag3D"},
                                                              {"fields_2d", "Fields/Diag2D"},
                                                              {"density_2d", "Particles/Electrons/Diag2D"},
                                                              {"energy_spectrum", "Particles/Electrons"},
                                                              {"probes", "Probes"}};
    for (const auto& [type, obstruction] : cases) {
        check::scenario(type + ": an ordinary file cannot serve as an output directory", [&] {
            const auto work = root / ("blocked_directory_" + type);
            fs::create_directory(work);
            fs::current_path(work);
            if (!obstruction.parent_path().empty())
                fs::create_directories(obstruction.parent_path());
            {
                std::ofstream file(obstruction);
                file << "preserve this file";
            }
            require_config_failure({{"outputs", nlohmann::json::array({{{"type", type}, {"interval", 1}}})}},
                                   obstruction.filename().string(), domain, species, field);
            std::ifstream file(obstruction);
            std::string contents;
            std::getline(file, contents);
            check::require(contents == "preserve this file", "failed preparation overwrote the existing file");
            check::require(!fs::exists("Recovery/Fields/FieldE.backup"),
                           "failed preparation wrote a partial checkpoint");
        });
    }
    for (const std::string type : {"energy_balance", "boundary_stats", "probes"}) {
        check::scenario(type + ": failure to open the output file is reported", [&] {
            const auto work = root / ("blocked_file_" + type);
            fs::create_directory(work);
            fs::current_path(work);
            const auto path = type == "probes" ? fs::path("Probes/probes.txt")
                                               : fs::path(type == "energy_balance" ? "energy.txt" : "boundary.txt");
            fs::create_directories(path);
            require_config_failure({{"outputs", nlohmann::json::array({{{"type", type}}})}}, "Cannot open", domain,
                                   species, field);
            check::require(fs::is_directory(path), "failed open changed the obstructing directory");
        });
    }
}

void supported_edge_settings_remain_valid(const fs::path& root, const Domain& domain, const Species& species,
                                          const Field3d& field) {
    const auto work = root / "valid_edge_settings";
    fs::create_directory(work);
    fs::current_path(work);
    const nlohmann::json recovery = {{"outputs", nlohmann::json::array({{{"type", "recovery"}}})}};
    check::require(OutputFactory::create(recovery, domain, species, {{"Dt", 1}}, field, field, field, field).empty(),
                   "absent recovery interval must keep recovery disabled");
    check::require(fs::is_empty(work), "disabled recovery created a directory");
    const auto settings = nlohmann::json::parse(R"({"outputs":[
        {"type":"probes","points":[[0,0,0],[4,4,4]],"species":"all"},
        {"type":"fields_2d","planes":{"z":[-100,100]},"interval":1,"average_steps":2}
    ]})");
    Diagnostics measured(species);
    auto outputs = OutputFactory::create(settings, domain, species, {{"Dt", 1}}, field, field, field, field);
    BoundaryConditionHandler boundaries;
    write_diagnostics(0, measured, {domain, species, field, field, field, field, boundaries}, outputs);
    for (auto& output : outputs) output->finalize();
    std::ifstream probes("Probes/probes.txt");
    std::string line;
    int rows = 0;
    while (std::getline(probes, line)) ++rows;
    check::require(rows == 3, "probes on either physical wall must remain supported");
    check::require(fs::is_empty("Fields/Diag2D"), "out-of-range planes must retain the skip behavior");
}
}   // namespace

void verify_diagnostic_output_failures(const fs::path& root) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    Species species;
    species.emplace("Electrons", std::make_unique<ParticlesArray>(check::species("Electrons"), domain));
    Field3d field(domain.size(), 3);
    field.setZero();
    check::scenario("All plane axes preserve coordinates and component ordering",
                    [&] { plane_axes_preserve_coordinates(root); });
    check::scenario("Buffered output failures are reported with filenames",
                    [&] { failed_writes_are_reported(root, domain, species, field); });
    invalid_settings_create_no_output(root, domain, species, field);
    blocked_directories_fail_before_writing(root, domain, species, field);
    check::scenario("Disabled recovery, wall probes and out-of-range planes retain their contracts",
                    [&] { supported_edge_settings_remain_valid(root, domain, species, field); });
}
