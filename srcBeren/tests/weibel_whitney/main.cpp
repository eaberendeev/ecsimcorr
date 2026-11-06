#include <bit>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>

#include "diagnostics/DensityMeasurements.h"
#include "simulation/ecsim/simulation_ecsim.h"
#include "simulation/ecsim_corr/simulation_ecsim_corr.h"
#include "utils/sparse.h"
#include "utils/containers.h"
#include "runtime/profiler.h"

namespace {
constexpr int modes[8][2] = {{1, 1}, {1, 0}, {0, 1}, {2, 0}, {0, 2}, {2, 1}, {1, 2}, {0, 3}};

// E and B share the code normalization: energy = 1/2 sum_active |field|^2.
double calc_energy_field(const Field3d& field, const IndexRange& range) {
    return 0.5 * dot_product_sum(field, field, range);
}

Field3d deposited_current(SimulationEcsim& s) {
    return dynamic_cast<SimulationEcsimCorr&>(s).fieldJe;
}

Field3d charge(SimulationEcsim& s) {
    return measure_charge_density(s.species, s.domain, s.bc_handler, ShapeType::Linear);
}

Field3d gauss(SimulationEcsim& s, const Field3d& rho) {
    const auto range = s.bc_handler.active_range(s.domain.grid);
    double mean = 0;
    for (int i = range.start[0]; i < range.end[0]; ++i)
        for (int j = range.start[1]; j < range.end[1]; ++j)
            for (int k = range.start[2]; k < range.end[2]; ++k) mean += rho(i, j, k, 0);
    mean /= s.domain.num_cells().elements_product();
    Field3d defect(s.domain.size(), 1);
    spmv(s.maxwell.divE, s.fieldEn, defect);
    // The single electron species has a uniform immobile neutralizing
    // background. Remove only its mean charge, not the CIC fluctuations.
    for (int i = 0; i < defect.capacity(); ++i) defect[i] -= rho[i] - mean;
    return defect;
}

void observe(SimulationEcsim& s, int step, int interval, const Field3d& rho, const Field3d& previous_charge,
             const Field3d& initial_defect, const Field3d& integrated_current, double initial_energy,
             std::ofstream& csv) {
    if (!s.fieldEn.data().allFinite() || !s.fieldBn.data().allFinite() || !rho.data().allFinite())
        throw std::runtime_error("non-finite Weibel fields or density");
    const auto range = s.bc_handler.active_range(s.domain.grid);
    const auto defect = gauss(s, rho);
    Field3d continuity(s.domain.size(), 1);
    spmv(s.maxwell.divE, integrated_current, continuity);
    continuity += rho - previous_charge;
    if (interval)
        continuity *= 1 / (interval * s.dt_);
    double gr2 = 0, gd2 = 0, dr2 = 0, cr2 = 0;
    std::complex<double> bx[8]{}, by[8]{};
    const double two_pi = 2 * std::acos(-1.0);
    const auto cells = s.domain.num_cells();
    const double count = cells.elements_product();
    for (int i = range.start[0]; i < range.end[0]; ++i)
        for (int j = range.start[1]; j < range.end[1]; ++j)
            for (int k = range.start[2]; k < range.end[2]; ++k) {
                gr2 += std::pow(defect(i, j, k, 0), 2);
                gd2 += std::pow(defect(i, j, k, 0) - initial_defect(i, j, k, 0), 2);
                dr2 += std::pow(rho(i, j, k, 0) + 1, 2);
                cr2 += std::pow(continuity(i, j, k, 0), 2);
                for (int m = 0; m < 8; ++m) {
                    const double phase = -two_pi * (modes[m][0] * double(i - range.start[0]) / cells[0] +
                                                    modes[m][1] * double(j - range.start[1]) / cells[1]);
                    const auto factor = std::polar(1.0 / count, phase);
                    bx[m] += factor * s.fieldBn(i, j, k, 0);
                    by[m] += factor * s.fieldBn(i, j, k, 1);
                }
            }
    const double we = calc_energy_field(s.fieldEn, range), wb = calc_energy_field(s.fieldBn, range);
    double kinetic = 0;
    for (const auto& [name, sp] : s.species) kinetic += sp->get_kinetic_energy();
    if (!std::isfinite(we + wb + kinetic + gr2 + gd2 + dr2 + cr2))
        throw std::runtime_error("non-finite Weibel observable");
    csv << step * s.dt_ << ',' << we << ',' << wb << ',' << kinetic << ','
        << (we + wb + kinetic - initial_energy) / initial_energy << ',' << std::sqrt(gr2 / count) << ','
        << std::sqrt(gd2 / count) << ',' << std::sqrt(dr2 / count) << ',' << std::sqrt(cr2 / count);
    for (int m = 0; m < 8; ++m) csv << ',' << std::sqrt(std::norm(bx[m]) + std::norm(by[m]));
    csv << ',' << std::arg(bx[0]) << '\n';
    csv.flush();
    if (!csv)
        throw std::runtime_error("cannot write Weibel observables");
}
}   // namespace

int main(int argc, char** argv) {
    try {
        nlohmann::json system, particles;
        std::ifstream("system_config.json") >> system;
        std::ifstream("particles_config.json") >> particles;
        if (system.value("StartFromTime", 0) != 0 || !system.at("Collisions").empty())
            throw std::invalid_argument("Weibel comparison requires fresh collisionless runs");
        auto simulation = build_simulation(system, particles, argc, argv);
        simulation->init();
        auto& s = dynamic_cast<SimulationEcsim&>(*simulation);
        // Order-independent fingerprint verifies identical individual particles
        // despite nondeterministic parallel insertion into cell buckets.
        uint64_t fingerprint = 0, population = 0;
        double initial_energy = 0;
        for (const auto& [name, sp] : s.species) {
            initial_energy += sp->get_kinetic_energy();
            for (int cell = 0; cell < sp->size(); ++cell)
                for (const auto& p : (*sp)(cell)) {
                    uint64_t hash = 1469598103934665603ULL;
                    for (int a = 0; a < 3; ++a) {
                        hash = (hash ^ std::bit_cast<uint64_t>(p.coord[a])) * 1099511628211ULL;
                        hash = (hash ^ std::bit_cast<uint64_t>(p.velocity[a])) * 1099511628211ULL;
                    }
                    fingerprint += hash;
                    ++population;
                }
        }
        std::ofstream("initial_state.json") << nlohmann::json{
            {"fingerprint", std::to_string(fingerprint)},
            {"particles", population},
            {"energy", initial_energy}}.dump(2);
        const int cadence = system.value("ComparisonInterval", 6), last = system.at("LastTimestep");
        const int profile_step = system.value("ComparisonProfileStep", 0);
        if (cadence <= 0 || last <= 0)
            throw std::invalid_argument("invalid Weibel comparison interval or duration");
        Field3d previous_charge = charge(s), initial_defect = gauss(s, previous_charge);
        Field3d integrated_current = Field3d::Zero(s.domain.size(), 3);
        std::ofstream csv("comparison.csv");
        csv << std::setprecision(17)
            << "time,electric,magnetic,kinetic,relative_energy_error,gauss_rms,gauss_drift_rms,density_rms,continuity_"
               "interval_rms";
        for (const auto& mode : modes) csv << ",mode" << mode[0] << mode[1];
        csv << ",phase11\n";
        s.make_diagnostic(0);
        observe(s, 0, 0, previous_charge, previous_charge, initial_defect, integrated_current, initial_energy, csv);
        const auto start = std::chrono::steady_clock::now();
        int previous_step = 0;
        for (int step = 1; step <= last; ++step) {
            s.prepare_step(step);
            s.make_step(step);
            // All comparison configurations have E_external=0. Corr retains
            // its full-step trajectory current in fieldJe; ordinary ECSIM
            // reconstructs the current used by the Maxwell solve.
            const std::string scheme = s.system_config.value("Scheme", std::string());
            if (scheme == "ecsim_corr")
                integrated_current += s.dt_ * deposited_current(s);
            else
                integrated_current += s.dt_ * s.fieldJp + 2.0 * (s.field_system.L * s.fieldEp);
            s.make_diagnostic(step);
            if (step % cadence == 0 || step == last) {
                const auto rho = charge(s);
                observe(s, step, step - previous_step, rho, previous_charge, initial_defect, integrated_current,
                        initial_energy, csv);
                previous_charge = rho;
                previous_step = step;
                integrated_current.setZero();
            }
            if (step == profile_step) {
                std::ofstream profile("step_profile.txt");
                timer::printTreeProfile(profile);
                timer::writeFullProfile("step_profile.json");
            }
            timer::clearFullProfile();
            timer::clearTree();
            if (step % 20 == 0 || step == last) {
                const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                std::cout << "progress " << step << '/' << last << " t=" << step * s.dt_ << " seconds=" << elapsed
                          << std::endl;
            }
        }
        s.finalize();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Weibel comparison failed: " << error.what() << '\n';
        return 1;
    }
}
