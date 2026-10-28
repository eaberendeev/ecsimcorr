#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "simulation.h"

int main(int argc, char** argv) {
    try {
        const int warmup = argc > 1 ? std::stoi(argv[1]) : 2;
        const int count = argc > 2 ? std::stoi(argv[2]) : 10;
        if (warmup < 0 || count < 1)
            throw std::invalid_argument("expected nonnegative warm-up and positive sample count");
        nlohmann::json system, particles;
        std::ifstream("system_config.json") >> system;
        std::ifstream("particles_config.json") >> particles;
        auto simulation = build_simulation(system, particles, argc, argv);
        simulation->init();
        simulation->make_diagnostic(0);
        std::ofstream csv("step_timing.csv");
        csv << "step,warmup,prepare_seconds,advance_seconds,diagnostic_seconds\n" << std::setprecision(12);
        std::vector<double> advance, total, diagnostic;
        using Clock = std::chrono::steady_clock;
        for (int step = 1; step <= warmup + count; ++step) {
            const auto start = Clock::now();
            simulation->prepare_step(step);
            const auto prepared = Clock::now();
            simulation->make_step(step);
            const auto advanced = Clock::now();
            simulation->make_diagnostic(step);
            const auto done = Clock::now();
            const double p = std::chrono::duration<double>(prepared - start).count();
            const double a = std::chrono::duration<double>(advanced - prepared).count();
            const double d = std::chrono::duration<double>(done - advanced).count();
            csv << step << ',' << (step <= warmup) << ',' << p << ',' << a << ',' << d << '\n';
            csv.flush();
            if (step > warmup) {
                advance.push_back(p + a);
                total.push_back(p + a + d);
                diagnostic.push_back(d);
            }
        }
        timer::writeFullProfile("step_profile.json");
        simulation->finalize();
        auto mean = [](const std::vector<double>& values) {
            return std::accumulate(values.begin(), values.end(), 0.) / values.size();
        };
        nlohmann::json report = {{"scheme", system.at("Scheme")},
                                 {"dt", system.at("Dt")},
                                 {"threads", omp_get_max_threads()},
                                 {"warmup_steps", warmup},
                                 {"measured_steps", count},
                                 {"mean_step_seconds", mean(advance)},
                                 {"mean_diagnostic_seconds", mean(diagnostic)},
                                 {"mean_total_seconds", mean(total)}};
        std::ofstream("step_timing.json") << report.dump(2) << '\n';
        std::cout << report.dump(2) << std::endl;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
