#include <omp.h>

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "Particle.h"
#include "collisions_with_neutrals.h"
#include "particles_distribution_collection.h"

// физические константы по умолчанию
const double DEFAULT_n0 = 1.e13;
const double DEFAULT_me = 1;
const double DEFAULT_mi = 1836;
const double DEFAULT_mn = DEFAULT_mi + 1;

// Параметры одного теста
struct TestCase {
    std::string name_prefix;
    double n0 = DEFAULT_n0;
    int num_particles = 1250000;
    long num_steps = 100000;
    double neutrals_energy_kev = 15.0;
    double particles_energy_kev = 1.0;
    double ncp = 0.1;
    double charged_mass = DEFAULT_me;
    double neutral_mass = DEFAULT_mn;
    uint64_t seed = 13;
    bool write_profile = true;
    double dt = 300;   // dt для каждого теста
    CollisionScheme scheme =
        CollisionScheme::PHYSICAL_ONLY;   // Значение по умолчанию
    CollisionProcessOptions process_opts =
        CollisionProcessOptions();   // Параметры столкновений
};

// Функция загрузки параметров теста из JSON
bool load_test_cases_from_json(const std::string& filename,
                               std::vector<TestCase>& test_cases) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open JSON file: " << filename << std::endl;
        return false;
    }

    nlohmann::json j;
    file >> j;

    if (j.find("tests") == j.end()) {
        std::cerr << "No 'tests' section found in " << filename << std::endl;
        return false;
    }

    for (const auto& test_data : j["tests"]) {
        TestCase tc;
        tc.name_prefix = test_data.value("name_prefix", tc.name_prefix);
        tc.num_particles = test_data.value("num_particles", tc.num_particles);
        tc.n0 = test_data.value("n0", tc.n0);
        tc.num_steps = test_data.value("num_steps", tc.num_steps);
        tc.neutrals_energy_kev =
            test_data.value("neutrals_energy_kev", tc.neutrals_energy_kev);
        tc.particles_energy_kev =
            test_data.value("particles_energy_kev", tc.particles_energy_kev);
        tc.ncp = test_data.value("ncp", tc.ncp);
        tc.charged_mass = test_data.value("charged_mass", tc.charged_mass);
        tc.neutral_mass = test_data.value("neutral_mass", tc.neutral_mass);
        tc.seed = test_data.value("seed", tc.seed);
        tc.write_profile = test_data.value("write_profile", tc.write_profile);
        tc.dt = test_data.value("dt", tc.dt);

        // Чтение CollisionScheme
        std::string scheme_str =
            test_data.value("collision_scheme", "PHYSICAL_ONLY");
        if (scheme_str == "PHYSICAL_ONLY") {
            tc.scheme = CollisionScheme::PHYSICAL_ONLY;
        } else if (scheme_str == "NULL_COLLISION") {
            tc.scheme = CollisionScheme::NULL_COLLISION;
        } else {
            std::cerr << "Unknown CollisionScheme: " << scheme_str << std::endl;
            return false;
        }

        // Чтение параметров столкновений
        auto collision_options = test_data["collision_options"];
        tc.process_opts.electron_ionization =
            collision_options.value("electron_ionization", false);
        tc.process_opts.proton_charge_exchange =
            collision_options.value("proton_charge_exchange", false);
        tc.process_opts.proton_ionization =
            collision_options.value("proton_ionization", false);

        test_cases.push_back(tc);
    }

    return true;
}

double velocity_from_kev(double kev, double mass) {
    return sqrt(2 * kev / 511.0 / mass);
}

// run_test: инициализация + основной цикл
bool run_test(const TestCase& tc) {
    const int num_particles = tc.num_particles;
    const long num_steps = tc.num_steps;
    const double neutrals_energy = tc.neutrals_energy_kev;
    const double particles_energy = tc.particles_energy_kev;
    const double ncp = tc.ncp;
    const double m_charged = tc.charged_mass;
    const double m_neutral = tc.neutral_mass;
    const double dt = tc.dt;

    // инициализация частиц
    double3 velocity_neutral = {velocity_from_kev(neutrals_energy, m_neutral),
                                0.0, 0.0};

    int initial_neutrals =
        static_cast<int>(std::max(0.0, double(num_particles) * ncp));
    std::vector<Particle> charged(num_particles);
    std::vector<Particle> neutrals(initial_neutrals);

    // RNG/распределения — инициализируем по seed
    ThreadRandomGenerator gen(static_cast<unsigned int>(tc.seed));

    // sigma для распределения скоростей заряженных частиц
    const double sigma = sqrt(particles_energy / 511.0 / m_charged);
    auto vel_dist = std::make_shared<GaussianVelocity>(double3(0.,0.,0), double3(sigma, sigma, sigma));
    for (auto & cp : charged){
        cp.velocity = vel_dist->sample(gen);
    }
    for(auto& np : neutrals){
        np.velocity = velocity_neutral;
    }

    ColliderWithNeutrals colliderWithNeutrals(tc.n0, tc.scheme,
                                              tc.process_opts);

    double freq_max = 0.0;

    std::string filename =
        tc.name_prefix + "_dt_" + std::to_string(dt) + ".txt";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }

    // Заголовок
    outfile << std::scientific << std::setprecision(6);
    outfile << "# step time " << "charged_count " <<  "charged_energy neutral_count neutral_energy\n";

    for (long step = 0; step < num_steps; ++step) {
        int current_neutral_count = static_cast<int>(neutrals.size());
        if(current_neutral_count < 0.05 * initial_neutrals) continue;

        std::uniform_int_distribution<int> dis(
            0, std::max(0, current_neutral_count - 1));

        for (int i = 0;
             i < static_cast<int>(charged.size()) && current_neutral_count > 0;
             ++i) {
            double3 vcp = charged[i].velocity;

            if (current_neutral_count != static_cast<int>(dis.max()) + 1) {
                dis = std::uniform_int_distribution<int>(
                    0, current_neutral_count - 1);
            }

            int randomIndex = dis(gen.gen());
            double3 vn = neutrals[randomIndex].velocity;
            double nn = double(current_neutral_count) / double(num_particles);

            auto [is_collided, ve_new, vi_new] =
                colliderWithNeutrals.collision_with_neutral(
                    vcp, vn, m_charged, m_neutral, 1.0, nn, dt, freq_max);

            bool remove_neutral = false;
            if (is_collided) {
                remove_neutral = is_collided;
            } else {
                // если vn изменился — значит collision обработал neutral
                // (например, charge-exchange)
                if (vn != neutrals[randomIndex].velocity) {
                    remove_neutral = true;
                }
            }

            if (remove_neutral) {
                std::swap(neutrals[randomIndex],
                          neutrals[current_neutral_count - 1]);
                --current_neutral_count;
            }
        }

        if (current_neutral_count > 0) {
            neutrals.resize(current_neutral_count);
        }

        // вычисление энергий
        double charged_energy = get_energy_particles(charged, m_charged, 1);
        double neutral_energy = get_energy_particles(neutrals, m_neutral, 1);

        double current_time = step * dt;
        outfile << step << " " << current_time << " " << charged.size() << " "
                << charged_energy << " " << neutrals.size() << " "
                << neutral_energy << "\n";

        // лог прогресса + профиль каждые 10 шагов. FEATURE: for parallel delete it or write to files
#pragma omp master
        if (step % 10 == 0) {
            std::cout << "Test [" << tc.name_prefix << "] Step " << step << "/"
                      << num_steps << ", time: " << current_time << ", "
                      << " charged count: " << charged.size() << ", energy "
                      << " charged energy: " << charged_energy
                      << ", neutrals: " << neutrals.size()
                      << " energy neutrals: " << neutral_energy << std::endl;
            std::cout << "=== Collision profiling summary for "
                      << tc.name_prefix << " ===\n";
            colliderWithNeutrals.profiler.print_report(std::cout);
            colliderWithNeutrals.profiler.reset();
        }
        
    }

    outfile.close();
    std::cout << "Data saved to: " << filename << std::endl;

    if (tc.write_profile) {
        std::string profname =
            "profile_" + tc.name_prefix + "_dt_" + std::to_string(dt) + ".txt";
        std::ofstream pf(profname);
        if (pf.is_open()) {
            colliderWithNeutrals.profiler.print_report(pf);
            pf.close();
            std::cout << "Profiler saved to: " << profname << std::endl;
        } else {
            std::cerr << "Cannot open profiler file: " << profname << std::endl;
        }
    }

    return true;
}

int main() {
    // Загрузка всех тестов из JSON
    std::string config_file = "test_config.json";
    std::vector<TestCase> test_cases;

    if (!load_test_cases_from_json(config_file, test_cases)) {
        return 1;
    }

    // Ограничиваем число потоков
    int max_threads = omp_get_num_threads();

    // Основная петля по тестам: создаём параллельные задачи OpenMP для каждого
    // теста
#pragma omp parallel num_threads(max_threads)
    {
#pragma omp single
        {
            for (const auto& tc : test_cases) {
#pragma omp task firstprivate(tc)
                { run_test(tc); }
            }
        }
    }

    return 0;
}
