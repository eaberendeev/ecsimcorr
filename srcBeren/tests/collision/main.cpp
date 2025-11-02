#include <omp.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "Particle.h"
#include "collisions_with_neutrals.h"
#include "particles_distribution_collection.h"

// физические константы по умолчанию (можно переопределять в тестах)
const double DEFAULT_n0 = 1.e13;
const double DEFAULT_me = 1;
const double DEFAULT_mi = 100;
const double DEFAULT_mn = DEFAULT_mi + 1;

double velocity_from_kev(double kev, double mass) {
    return sqrt(2 * kev / 511.0 / mass);
}

// Параметры одного теста
struct TestCase {
    std::string name_prefix;   // префикс для файла и логов
    int num_particles = 1000000;
    long num_steps = 100;
    double neutrals_energy_kev = 15.0;
    double particles_energy_kev = 1.0;   // electrons or ions energy
    double ncp = 0.1;   // neutrals per charged particle (fraction)
    double charged_mass = DEFAULT_me;
    double neutral_mass = DEFAULT_mn;
    CollisionScheme scheme = CollisionScheme::PHYSICAL_ONLY;
    CollisionProcessOptions process_opts = CollisionProcessOptions();
    bool use_is_collided = true;   // если true — используем is_collided to
                                   // remove neutral, иначе compare vn change
    std::string particle_label = "particles";
    uint64_t seed = 13;   // seed для воспроизводимости
    bool write_profile =
        false;   // если true, сохранить профиль в файл по окончании прогона
    TestCase() = default;
};

// run_test: инициализация + основной цикл (как в старой непакетной версии).
// Важно: все ресурсы внутри run_test локальны (каждый тест создаёт свой
// ColliderWithNeutrals, свой ThreadRandomGenerator и свои векторы частиц),
// поэтому безопасно вызывать из OpenMP-задачи.
bool run_test(const TestCase& tc, double dt) {
    const int num_particles = tc.num_particles;
    const long num_steps = tc.num_steps;
    const double neutrals_energy = tc.neutrals_energy_kev;
    const double particles_energy = tc.particles_energy_kev;
    const double ncp = tc.ncp;
    const double m_charged = tc.charged_mass;
    const double m_neutral = tc.neutral_mass;

    // инициализация частиц
    double3 velocity_neutral = {velocity_from_kev(neutrals_energy, m_neutral),
                                0.0, 0.0};

    int initial_neutrals =
        static_cast<int>(std::max(0.0, double(num_particles) * ncp));
    std::vector<Particle> charged(num_particles);
    std::vector<Particle> neutrals(initial_neutrals);

    // RNG/распределения — инициализируем по seed, чтобы тест был воспроизводим
    ThreadRandomGenerator gen(static_cast<unsigned int>(tc.seed));

    // sigma для распределения скоростей заряженных частиц
    const double sigma = sqrt(particles_energy / 511.0 / m_charged);
    distribute_pulse_gauss(charged, double3(sigma, sigma, sigma), gen);
    set_velocity(neutrals, velocity_neutral);

    ColliderWithNeutrals colliderWithNeutrals(DEFAULT_n0, tc.scheme,
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
    outfile << "# step time " << tc.particle_label << "_count "
            << tc.particle_label << "_energy neutral_count neutral_energy\n";

    for (long step = 0; step < num_steps; ++step) {
        int current_neutral_count = static_cast<int>(neutrals.size());

        // создаём распределение для индекса нейтрала — один раз, и обновляем
        // только при изменении размера
        std::uniform_int_distribution<int> dis(
            0, std::max(0, current_neutral_count - 1));

        for (int i = 0;
             i < static_cast<int>(charged.size()) && current_neutral_count > 0;
             ++i) {
            double3 vcp = charged[i].velocity;

            // если число нейтралов изменилось, обновляем distribution bound
            if (current_neutral_count != static_cast<int>(dis.max()) + 1) {
                dis = std::uniform_int_distribution<int>(
                    0, current_neutral_count - 1);
            }

            int randomIndex = dis(gen.gen());   // gen.gen() -> быстрый uint64_t
                                                // source
            double3 vn = neutrals[randomIndex].velocity;
            double nn = double(current_neutral_count) / double(num_particles);

            auto [is_collided, ve_new, vi_new] =
                colliderWithNeutrals.collision_with_neutral(
                    vcp, vn, m_charged, m_neutral, 1.0, nn, dt, freq_max);

            bool remove_neutral = false;
            if (tc.use_is_collided) {
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

            // В оригинале charged[i] не обновлялся — сохраняем это поведение.
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
        if (step % 10 == 0) {
            std::cout << "Test [" << tc.name_prefix << "] Step " << step << "/"
                      << num_steps << ", time: " << current_time << ", "
                      << tc.particle_label << ": " << charged.size()
                      << ", energy " << tc.particle_label << ": "
                      << charged_energy << ", neutrals: " << neutrals.size()
                      << " energy neutrals: " << neutral_energy << std::endl;

            std::cout << "=== Collision profiling summary for "
                      << tc.name_prefix << " ===\n";
            colliderWithNeutrals.profiler.print_report(std::cout);
            colliderWithNeutrals.profiler.reset();
        }
    }   // end steps

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
    // подготовка трёх тестов: electron ionization, proton ionization, charge
    // exchange

    TestCase tc_electron;
    tc_electron.name_prefix = "test_electron_ionization";
    tc_electron.particle_label = "electrons";
    tc_electron.charged_mass = DEFAULT_me;
    tc_electron.neutral_mass = DEFAULT_mn;
    tc_electron.particles_energy_kev = 1.0;
    tc_electron.neutrals_energy_kev = 15.0;
    tc_electron.ncp = 0.1;
    tc_electron.scheme = CollisionScheme::PHYSICAL_ONLY;
    tc_electron.process_opts = CollisionProcessOptions();
    tc_electron.process_opts.electron_ionization = true;
    tc_electron.process_opts.proton_charge_exchange = false;
    tc_electron.process_opts.proton_ionization = false;
    tc_electron.use_is_collided = true;
    tc_electron.seed = 42;
    tc_electron.write_profile = true;

    TestCase tc_proton;
    tc_proton.name_prefix = "test_proton_ionization";
    tc_proton.particle_label = "ions";
    tc_proton.charged_mass = DEFAULT_mi;
    tc_proton.neutral_mass = DEFAULT_mn;
    tc_proton.particles_energy_kev = 1.0;
    tc_proton.neutrals_energy_kev = 15.0;
    tc_proton.ncp = 0.1;
    tc_proton.scheme = CollisionScheme::PHYSICAL_ONLY;
    tc_proton.process_opts = CollisionProcessOptions();
    tc_proton.process_opts.electron_ionization = false;
    tc_proton.process_opts.proton_charge_exchange = false;
    tc_proton.process_opts.proton_ionization = true;
    tc_proton.use_is_collided = true;
    tc_proton.seed = 42;
    tc_proton.write_profile = true;

    TestCase tc_cx;
    tc_cx.name_prefix = "test_charge_exchange";
    tc_cx.particle_label = "ions";
    tc_cx.charged_mass = DEFAULT_mi;
    tc_cx.neutral_mass = DEFAULT_mn;
    tc_cx.particles_energy_kev = 1.0;
    tc_cx.neutrals_energy_kev = 15.0;
    tc_cx.ncp = 0.1;
    tc_cx.scheme = CollisionScheme::PHYSICAL_ONLY;
    tc_cx.process_opts = CollisionProcessOptions();
    tc_cx.process_opts.electron_ionization = false;
    tc_cx.process_opts.proton_charge_exchange = true;
    tc_cx.process_opts.proton_ionization = false;
    tc_cx.use_is_collided = false;   // сравнение скоростей нейтрала
    tc_cx.seed = 42;
    tc_cx.write_profile = true;

    std::vector<double> dt_values = {150, 15, 1.5};  // in 1/w_pe

    // Основная петля по dt: для каждого dt создаём три OpenMP tasks (по одному
    // на тест). Порядок: ждём завершения задач текущего dt (taskwait), затем
    // переходим к следующему dt.
#pragma omp parallel num_threads(1)
    {
#pragma omp single
        {
            for (double dt : dt_values) {
                std::cout << "Launching tests for dt = " << dt << std::endl;

                // task for electron
#pragma omp task firstprivate(dt, tc_electron)
                {
                    run_test(tc_electron, dt);
                }

                // task for proton
#pragma omp task firstprivate(dt, tc_proton)
                {
                    run_test(tc_proton, dt);
                }

                // task for charge-exchange
#pragma omp task firstprivate(dt, tc_cx)
                {
                    run_test(tc_cx, dt);
                }

                // wait for the three tasks to finish before starting next dt
                // iteration
#pragma omp taskwait
            }   // end for dt
        }   // end single
    }   // end parallel

    return 0;
}
