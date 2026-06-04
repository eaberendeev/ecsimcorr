#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "World.h"
#include "boundary_conditions.h"
#include "nlohmann/json.hpp"
#include "particles/ParticlesArray.h"
#include "simulation/ecsim_corr/algorithms_ecsim_corr.h"
#include "timer.h"

using Clock = std::chrono::high_resolution_clock;

inline double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct TimingAccum {
    double total_s = 0;
    int calls = 0;
    void add(double t) {
        total_s += t;
        ++calls;
    }
};

struct TestTimings {
    std::string name;
    int steps;
    double dt;
    int init_particles;
    std::unordered_map<std::string, TimingAccum> funcs;
};

std::vector<TestTimings> all_timings;

Domain create_domain_from_config(const nlohmann::json& cfg) {
    Domain dom;
    const nlohmann::json& d = cfg["domain"];
    const nlohmann::json& cs = d["cell_size"];
    const nlohmann::json& nc = d["num_cells"];
    dom.init(Vector3I(nc[0].get<int>(), nc[1].get<int>(), nc[2].get<int>()),
             Vector3R(cs[0].get<double>(), cs[1].get<double>(), cs[2].get<double>()));
    if (d.contains("cylinder")) {
        const nlohmann::json& cyl = d["cylinder"];
        dom.set_cylinder(Vector3R(cyl["center"][0].get<double>(), cyl["center"][1].get<double>(), 0),
                         cyl["radius"].get<double>());
    }
    return dom;
}

BoundaryConditionHandler create_bc_handler(const nlohmann::json& cfg, const Domain& domain) {
    BoundaryConditionHandler handler;
    if (cfg.contains("boundary_conditions")) {
        nlohmann::json bc_cfg;
        bc_cfg["Boundary_conditions"] = cfg["boundary_conditions"];
        handler.load_from_json(bc_cfg, domain);
    }
    return handler;
}

bool run_multi_step_test(const nlohmann::json& config) {
    const std::string name = config["name"];
    std::cout << "Testing " << name << std::endl;

    Domain domain = create_domain_from_config(config);
    BoundaryConditionHandler bc_handler = create_bc_handler(config, domain);

    const int steps = config.value("steps", 1);
    const double dt = config.value("dt", 0.1);
    const double tolerance = config.value("tolerance", 1e-12);

    ParticlesArray ref_arr(config["species"], domain);
    ref_arr.distribute_initial_particles(ref_arr.get_initial_distributions(), domain);

    ParticlesArray test_arr(config["species"], domain);
    test_arr.distribute_initial_particles(test_arr.get_initial_distributions(), domain);

    TestTimings timings;
    timings.name = name;
    timings.steps = steps;
    timings.dt = dt;
    timings.init_particles = ref_arr.get_total_num_of_particles();

    std::cout << "  Particles: " << timings.init_particles << ", steps: " << steps << ", dt: " << dt << std::endl;

    Species empty_species;

    Field3d J_ref(domain.size(), 3);
    Field3d J_test(domain.size(), 3);

    bool all_ok = true;

    for (int step = 0; step < steps; ++step) {
        J_ref.setZero();
        J_test.setZero();

        auto t0 = Clock::now();
        move_and_calc_current_reference(ref_arr, dt, J_ref);
        timings.funcs["move_and_calc_current_ref"].add(elapsed(t0));

        t0 = Clock::now();
        move_and_calc_current(test_arr, dt, J_test);
        timings.funcs["move_and_calc_current"].add(elapsed(t0));

        const double errJ = (J_ref - J_test).norm();
        const double normJ = J_ref.norm();
        const double rel_errJ = normJ > 0 ? errJ / normJ : errJ;
        std::cout << "  Step " << step << ": ||J_ref - J_test|| = " << errJ;
        if (normJ > 0)
            std::cout << " (rel: " << rel_errJ << ")";
        std::cout << std::endl;

        if (rel_errJ >= tolerance) {
            std::cerr << "  FAIL: J mismatch at step " << step << std::endl;
            all_ok = false;
        }

        t0 = Clock::now();
        ref_arr.update_cells(domain);
        timings.funcs["update_cells_ref"].add(elapsed(t0));

        t0 = Clock::now();
        test_arr.update_cells(domain);
        timings.funcs["update_cells_test"].add(elapsed(t0));

        t0 = Clock::now();
        bc_handler.apply_to_particles(ref_arr, empty_species, domain);
        timings.funcs["apply_to_particles_ref"].add(elapsed(t0));

        t0 = Clock::now();
        bc_handler.apply_to_particles(test_arr, empty_species, domain);
        timings.funcs["apply_to_particles_test"].add(elapsed(t0));

        int count_ref = ref_arr.get_total_num_of_particles();
        int count_test = test_arr.get_total_num_of_particles();
        std::cout << "  Step " << step << " particles: ref=" << count_ref << " test=" << count_test << std::endl;

        // NOTE: for rare cases, this inequality could be hold for exactly the same function, since our parallel
        // functions does not fullfill run to run reproducibility
        if (count_ref != count_test) {
            std::cerr << "  FAIL: particle count mismatch at step " << step << std::endl;
            all_ok = false;
        }
    }

    if (all_ok) {
        auto t0 = Clock::now();
        ref_arr.density_on_grid_update_reference();
        timings.funcs["density_on_grid_update_ref"].add(elapsed(t0));

        t0 = Clock::now();
        test_arr.density_on_grid_update();
        timings.funcs["density_on_grid_update"].add(elapsed(t0));

        const double errDens = (ref_arr.densityOnGrid - test_arr.densityOnGrid).norm();
        const double normDens = ref_arr.densityOnGrid.norm();
        const double rel_errDens = normDens > 0 ? errDens / normDens : errDens;
        std::cout << "  Density: ||ref - test|| = " << errDens;
        if (normDens > 0)
            std::cout << " (rel: " << rel_errDens << ")";
        std::cout << std::endl;

        if (rel_errDens >= tolerance) {
            std::cerr << "  FAIL: density mismatch" << std::endl;
            all_ok = false;
        }
    }

    all_timings.push_back(timings);
    return all_ok;
}

void write_timing_summary(const std::string& pathStr) {
    const std::filesystem::path path = std::filesystem::absolute(pathStr);
    nlohmann::json out;
    out["tests"] = nlohmann::json::array();

    for (const TestTimings& t : all_timings) {
        nlohmann::json entry;
        entry["name"] = t.name;
        entry["steps"] = t.steps;
        entry["dt"] = t.dt;
        entry["initial_particles"] = t.init_particles;
        entry["functions"] = nlohmann::json::object();

        for (const auto& [fname, acc] : t.funcs) {
            nlohmann::json finfo;
            finfo["total_s"] = acc.total_s;
            finfo["calls"] = acc.calls;
            finfo["avg_ms"] = acc.total_s * 1000.0 / acc.calls;
            entry["functions"][fname] = finfo;
        }

        out["tests"].push_back(entry);
    }

    std::ofstream fout(path);
    fout << out.dump(2) << std::endl;
    std::cout << "Timing summary written to " << path << std::endl;
}

int main() {
    std::filesystem::path inputFile = std::filesystem::absolute("config.json");
    nlohmann::json config;
    {
        std::ifstream fin(inputFile);
        if (!fin) {
            std::cerr << "failed to read input config file: " << inputFile << std::endl;
            return 1;
        }
        fin >> config;
    }

    bool result = true;

    for (const nlohmann::json& it : config["cases"]) {
        const std::string name = it["name"];
        bool res = run_multi_step_test(it);

        timer::printSlice(std::cout, "void move_and_calc_current_reference", "void move_and_calc_current",
                          "void ParticlesArray::update_cells", "void BoundaryConditionHandler::apply_to_particles",
                          "void ParticlesArray::density_on_grid_update_reference",
                          "void ParticlesArray::density_on_grid_update");

        timer::clearTree();

        std::cout << "Test " << name << ": " << (res ? "PASS" : "FAIL") << std::endl << std::endl;
        result = result && res;
    }

    write_timing_summary("particles_array_timing.json");
    timer::writeTimerTree("particles_array_profile.json");

    return result ? 0 : 1;
}
