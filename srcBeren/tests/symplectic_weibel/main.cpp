#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>

#include "simulation_symplectic.h"

// Quiet deterministic 1D3V plasma. The benchmark subclass only sets initial
// data; every time step is the production SimulationSymplectic implementation.
class WeibelSimulation : public SimulationSymplectic {
   public:
    WeibelSimulation(const nlohmann::json& config, const nlohmann::json& particles, const nlohmann::json& experiment,
                     int argc, char** argv)
        : SimulationSymplectic(config, particles, argc, argv), experiment_(experiment) {
    }

    void init_fields() override {
        SimulationSymplectic::init_fields();
        const double k = experiment_.at("k"), amplitude = experiment_.at("seed_b");
        const double gamma = experiment_.at("seed_gamma");
        const double discrete_k = 2 * std::sin(0.5 * k * domain.cell_size(0)) / domain.cell_size(0);
        const int g = domain.grid.ghost_cells();
        for (int i = 0; i < domain.num_cells(0); ++i) {
            const auto xb = domain.get_node_position(i + g, g, g, FieldType::MAGNETIC, 2);
            const auto xe = domain.get_node_position(i + g, g, g, FieldType::ELECTRIC, 1);
            magnetic(i + g, g, g, 2) = amplitude * std::cos(k * xb[0]);
            electric(i + g, g, g, 1) = -gamma / discrete_k * amplitude * std::sin(k * xe[0]);
        }
        symplectic::WhitneyPIC(domain).sync(electric);
        symplectic::WhitneyPIC(domain).sync(magnetic);
    }
    void init_particles(const nlohmann::json& config) override {
        Simulation::init_particles(config);
        auto& electrons = *species.at("Electrons");
        auto& ions = *species.at("Ions");
        const auto velocities = experiment_.at("velocities").get<std::vector<std::vector<double>>>();
        for (int i = 0; i < domain.num_cells(0); ++i) {
            const Vector3R x((i + 0.5) * domain.cell_size(0), 0.5 * domain.cell_size(1), 0.5 * domain.cell_size(2));
            for (const auto& v : velocities) electrons.add_particle(Particle(x, {v.at(0), v.at(1), v.at(2)}));
            // m_i/m_e=1e12 approximates immobile ions, while preserving exact
            // initial nodal neutrality without an implicit background charge.
            ions.add_particle(Particle(x, {0, 0, 0}));
        }
        for (auto& [name, sp] : species) {
            deposit_density(*sp);
            sp->kineticEnergy = sp->get_kinetic_energy();
        }
    }

   private:
    nlohmann::json experiment_;
};

int main(int argc, char** argv) {
    try {
        std::ifstream input("weibel.json");
        if (!input)
            throw std::runtime_error("missing weibel.json; use Weibel_symplectic/run.py");
        nlohmann::json experiment;
        input >> experiment;
        auto system = experiment.at("system"), particles = experiment.at("particles");
        WeibelSimulation sim(system, particles, experiment, argc, argv);
        sim.init();
        symplectic::WhitneyPIC pic(sim.domain);
        std::ofstream csv("weibel.csv");
        csv << "step,time,magnetic_energy,electric_energy,kinetic_energy,relative_energy_error,gauss_error,div_b,"
               "mode_b_real,mode_b_imag,mode_e_real,mode_e_imag,vx2,vy2,vz2,wall_seconds\n";
        csv << std::setprecision(17);
        const int nx = sim.domain.num_cells(0), g = sim.domain.grid.ghost_cells();
        const double k = experiment.at("k"), dt = system.at("Dt");
        const int last = system.at("LastTimestep"), interval = experiment.at("sample_interval");
        double initial_energy = 0;
        const auto started = std::chrono::steady_clock::now();
        for (int step = 0; step <= last; ++step) {
            if (step > 0) {
                sim.prepare_step(step);
                sim.make_step(step);
            }
            if (step % interval != 0 && step != last)
                continue;
            sim.make_diagnostic(step);
            const double wb = pic.field_energy(sim.magnetic), we = pic.field_energy(sim.electric);
            double kinetic = 0;
            for (const auto& [name, sp] : sim.species) kinetic += sp->get_kinetic_energy();
            const double energy = kinetic + we + wb;
            if (step == 0)
                initial_energy = energy;
            const double relative_error = (energy - initial_energy) / initial_energy;
            std::complex<double> bmode(0, 0), emode(0, 0);
            double gauss = 0, div_b = 0;
            Field3d rho(sim.domain.size(), 1);
            sim.collect_charge_density(rho);
            for (int i = 0; i < nx; ++i) {
                const auto xb = sim.domain.get_node_position(i + g, g, g, FieldType::MAGNETIC, 2);
                const auto xe = sim.domain.get_node_position(i + g, g, g, FieldType::ELECTRIC, 1);
                bmode += sim.magnetic(i + g, g, g, 2) * std::polar(2.0 / nx, -k * xb[0]);
                emode += sim.electric(i + g, g, g, 1) * std::polar(2.0 / nx, -k * xe[0]);
                gauss = std::max(gauss, std::abs(pic.divergence(sim.electric, i, 0, 0, false) - rho(i + g, g, g, 0)));
                div_b = std::max(div_b, std::abs(pic.divergence(sim.magnetic, i, 0, 0, true)));
            }
            Vector3R second_moment(0, 0, 0);
            const auto& electrons = *sim.species.at("Electrons");
            for (int cell = 0; cell < electrons.size(); ++cell)
                for (const auto& p : electrons.particlesData(cell))
                    for (int d = 0; d < 3; ++d) second_moment[d] += p.velocity[d] * p.velocity[d];
            second_moment /= static_cast<double>(electrons.get_total_num_of_particles());
            if (!std::isfinite(energy) || !std::isfinite(std::abs(bmode)) || !std::isfinite(gauss) ||
                !std::isfinite(div_b) || std::abs(relative_error) > 0.05 || gauss > 1e-9 || div_b > 1e-11)
                throw std::runtime_error(
                    "Weibel monitor rejected non-finite state, energy drift or constraint violation");
            const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            csv << step << ',' << step * dt << ',' << wb << ',' << we << ',' << kinetic << ',' << relative_error << ','
                << gauss << ',' << div_b << ',' << bmode.real() << ',' << bmode.imag() << ',' << emode.real() << ','
                << emode.imag() << ',' << second_moment[0] << ',' << second_moment[1] << ',' << second_moment[2] << ','
                << elapsed << '\n';
            csv.flush();
            if (step == 0 || step % (10 * interval) == 0 || step == last)
                std::cout << "step=" << step << '/' << last << " t=" << step * dt << " |B_k|=" << std::abs(bmode)
                          << " dH/H=" << relative_error << " Gauss=" << gauss << " wall=" << elapsed << std::endl;
        }
        sim.finalize();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
