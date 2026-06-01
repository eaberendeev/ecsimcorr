#include <chrono>
#include <filesystem>
#include <fstream>
#include <random>

#include "nlohmann/json.hpp"
#include "particles/ParticlesArray.h"
#include "simulation/ecsim_corr/algorithms_ecsim_corr.h"
#include "timer.h"

ParticlesArray createParticlesArray(const nlohmann::json& config, const Domain& domain) {
    ParticlesArray array(config, domain);

    const double init_energy = array.distribute_initial_particles(array.get_initial_distributions(), domain);
    std::cout << array.name() << " distibuted with init energy: " << init_energy << "\n";

    return array;
}

std::chrono::high_resolution_clock::time_point now() {
    return std::chrono::high_resolution_clock::now();
}

double dur(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end) {
    return std::chrono::duration<double>(end - start).count();
}

bool runTest(const nlohmann::json& config) {
    Domain domain;
    domain.init({60, 60, 40}, {0.5, 0.5, 0.5});

    ParticlesArray refArray = createParticlesArray(config, domain);
    ParticlesArray array = createParticlesArray(config, domain);

    Field3d refJ(domain.size(), 3);
    refJ.setZero();
    std::chrono::time_point timePoint1 = now();
    refArray.density_on_grid_update_reference();
    std::chrono::time_point timePoint2 = now();
    move_and_calc_current_reference(refArray, 1.5, refJ);
    std::chrono::time_point timePoint3 = now();

    const Field3d refDensity = refArray.densityOnGrid;

    Field3d J(domain.size(), 3);
    J.setZero();
    std::chrono::time_point timePoint4 = now();
    array.density_on_grid_update();
    std::chrono::time_point timePoint5 = now();
    move_and_calc_current(array, 1.5, J);
    std::chrono::time_point timePoint6 = now();
    const Field3d density = array.densityOnGrid;

    const double errJ = (refJ - J).norm();
    const double errDensity = (refDensity - density).norm();
    std::cout << "density_on_grid_update: reference: " << dur(timePoint1, timePoint2)
              << " test: " << dur(timePoint4, timePoint5) << "s" << std::endl;
    std::cout << "move_and_calc_current: reference: " << dur(timePoint2, timePoint3)
              << " test: " << dur(timePoint5, timePoint6) << "s" << std::endl;
    std::cout << "||J_ref - J||: " << errJ << " ||J_ref - J||/||j_ref||: " << errJ / refJ.norm() << std::endl;
    std::cout << "||density_ref - density||: " << errDensity
              << " ||density_ref - density||/||density_ref||: " << errDensity / refJ.norm() << std::endl;

    return std::max(errJ, errDensity) < 1e-14;
}

int main() {
    std::filesystem::path inputFile = std::filesystem::absolute("config.json");
    nlohmann::json config;
    std::ifstream fin(inputFile);
    if (!fin) {
        std::cerr << "failed to read input config file: " << inputFile << std::endl;
        return 1;
    }
    fin >> config;
    fin.close();

    bool result = true;

    for (const auto& it : config["cases"]) {
        const std::string name = it["name"];
        const nlohmann::json& testConfig = it["config"];
        std::cout << "Testing " << name << std::endl;
        const bool res = runTest(testConfig);
        std::cout << "Test " << name << ": " << (res ? "pass" : "fail") << std::endl;
        result = res && result;
    }

    timer::writeTimerTree("particles_array_profile.json");

    return result ? 0 : 1;
}
