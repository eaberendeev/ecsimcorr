#include <omp.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <string>
#include <vector>

#include "ParticlesArray.h"
#include "collision_manager.h"
#include "config/domainconfig.h"
#include "grid/domain.h"
#include "sgs.h"

// ============================================================================
// Statistics
// ============================================================================
// Temperatures are returned in energy units (k_B = 1) following Takizuka & Abe
// (J. Comput. Phys. 25, 205, 1977):
//   T_parallel = m <(vx-ux)^2>
//   T_perp     = 0.5 * m * (<(vy-uy)^2> + <(vz-uz)^2>)
//   T_total    = (T_parallel + 2 T_perp) / 3
// Velocities are in units of c, masses in units of m_e, so T is in units of
// m_e c^2 = 511 keV.
struct Stats {
    double T_par;
    double T_perp;
    double T_total;
    double Etotal;
};

Stats compute_stats(ParticlesArray& sp) {
    Vector3R mean_velocity(0.0);
    size_t count = 0;
    for (int cell = 0; cell < sp.size(); ++cell)
        for (const auto& particle : sp.particlesData(cell)) {
            mean_velocity += particle.velocity;
            ++count;
        }
    if (count == 0)
        return {0, 0, 0, 0};
    mean_velocity /= static_cast<double>(count);

    Vector3R thermal_second_moment(0.0);
    double kinetic_second_moment = 0;
    for (int cell = 0; cell < sp.size(); ++cell)
        for (const auto& particle : sp.particlesData(cell)) {
            const Vector3R thermal_velocity = particle.velocity - mean_velocity;
            for (int axis = 0; axis < 3; ++axis)
                thermal_second_moment[axis] += thermal_velocity[axis] * thermal_velocity[axis];
            kinetic_second_moment += particle.velocity.squaredNorm();
        }
    const double mass = sp.mass();
    const double parallel = mass * thermal_second_moment.x() / count;
    const double perpendicular = 0.5 * mass * (thermal_second_moment.y() + thermal_second_moment.z()) / count;
    // Temperature excludes the instantaneous bulk flow; conserved energy includes it.
    return {parallel, perpendicular, (parallel + 2 * perpendicular) / 3, 0.5 * sp.mpw() * mass * kinetic_second_moment};
}

void fill_maxwellian(ParticlesArray& sp, int num_particles, double sigma_x, double sigma_y, double sigma_z,
                     unsigned seed) {
    std::mt19937 rndEng(seed);
    std::normal_distribution<double> dist_x(0.0, sigma_x);
    std::normal_distribution<double> dist_y(0.0, sigma_y);
    std::normal_distribution<double> dist_z(0.0, sigma_z);

    // Full 3D Cartesian Maxwellian (vx, vy, vz independent Gaussians): sigma_x
    // sets T_parallel, sigma_y = sigma_z set T_perp for the isotropization test;
    // the equilibration tests pass the three sigmas equal (isotropic).  Subtract
    // the mean velocity so the net bulk flow is exactly zero (a spurious bulk
    // velocity would bias the measured temperature).
    std::vector<Vector3R> vels(num_particles);
    Vector3R mean(0.0, 0.0, 0.0);
    for (int i = 0; i < num_particles; i++) {
        vels[i] = Vector3R(dist_x(rndEng), dist_y(rndEng), dist_z(rndEng));
        mean += vels[i];
    }
    mean /= static_cast<double>(num_particles);
    for (int i = 0; i < num_particles; i++) {
        Particle p;
        p.coord = Vector3R(0.5, 0.5, 0.5);
        p.velocity = vels[i] - mean;
        sp.add_particle(p);
    }
}

nlohmann::json make_domain_config() {
    nlohmann::json config;
    config["Dx"] = 1.0;
    config["Dy"] = 1.0;
    config["Dz"] = 1.0;
    config["NumCellsX"] = 1;
    config["NumCellsY"] = 1;
    config["NumCellsZ"] = 1;
    return config;
}

nlohmann::json make_species_config(const std::string& name, double charge, double mass, int npc) {
    nlohmann::json config;
    config["Name"] = name;
    config["Charge"] = charge;
    config["Mass"] = mass;
    config["Density"] = 1.0;
    config["NumPartPerCell"] = npc;
    config["distribution"] = nlohmann::json::array();
    return config;
}

// Per-step variance <delta^2> at the thermal relative speed, as a sanity check:
// the Takizuka-Abe small-angle approximation requires <delta^2> << 1.
double estimate_delta2_per_step(double charge, double mass, double n_density, double u_thermal, double n0,
                                double coulomb_log, double dt) {
    const double wp3 = pow(SGS::get_plasma_freq(n0), 3);
    const double c3 = SGS::c * SGS::c * SGS::c;
    const double m_red = mass / 2.0;
    const double q = charge;
    return (wp3 / (c3 * n0)) * (coulomb_log * q * q * q * q * n_density * dt) /
           (8.0 * M_PI * m_red * m_red * u_thermal * u_thermal * u_thermal);
}

// ============================================================================
// Test 1: single-species temperature isotropization (anisotropic -> isotropic)
// ============================================================================
void run_test_isotropization(const std::string& test_name, double charge, double mass, double T_par0, double T_perp0,
                             int num_particles, int npc, int num_steps, double dt, double n0, double coulomb_log,
                             int diag_interval) {
    std::cout << "\n=== " << test_name << " (isotropization) ===\n";

    auto domain_config = make_domain_config();
    Domain domain = domain_from_config(domain_config);

    auto sp_config = make_species_config("TestSpecies", charge, mass, npc);
    ParticlesArray sp(sp_config, domain);

    // T_par = m sigma_x^2, T_perp = m sigma_yz^2  =>  sigma = sqrt(T/m).
    const double sigma_x = std::sqrt(T_par0 / mass);
    const double sigma_yz = std::sqrt(T_perp0 / mass);
    fill_maxwellian(sp, num_particles, sigma_x, sigma_yz, sigma_yz, 12345);

    CoulombCollisionOperator collider("TestSpecies", "TestSpecies", n0, coulomb_log);

    Stats s0 = compute_stats(sp);
    const double n_density = (double) num_particles / npc;
    const double u_th = std::sqrt(2.0 * s0.T_total / mass);   // thermal relative speed
    const double d2 = estimate_delta2_per_step(charge, mass, n_density, u_th, n0, coulomb_log, dt);
    std::cout << "  N=" << num_particles << " npc=" << npc << " n_density=" << n_density << " dt=" << dt
              << " steps=" << num_steps << "\n";
    std::cout << "  Initial: T_par=" << s0.T_par << " T_perp=" << s0.T_perp << " ratio=" << s0.T_perp / s0.T_par
              << " E=" << s0.Etotal << "\n";
    std::cout << "  <delta^2>/step (thermal) ~ " << d2 << " (should be << 1)\n";

    std::string fname = test_name + ".csv";
    std::ofstream out(fname);
    if (!out)
        throw std::runtime_error("Cannot write relaxation data: " + fname);
    out << std::setprecision(17);
    out << "step,time,T_par,T_perp,T_total,Etotal,Etotal_rel_err\n";
    out << "0,0," << s0.T_par << "," << s0.T_perp << "," << s0.T_total << "," << s0.Etotal << ",0\n";

    for (int step = 1; step <= num_steps; step++) {
        collider.collide_same_type(sp, dt);
        if (step % diag_interval == 0) {
            Stats s = compute_stats(sp);
            double rel_err = (s.Etotal - s0.Etotal) / s0.Etotal;
            out << step << "," << step * dt << "," << s.T_par << "," << s.T_perp << "," << s.T_total << "," << s.Etotal
                << "," << rel_err << "\n";
        }
    }
    out.close();

    nlohmann::json params;
    params["test_type"] = "isotropization";
    params["n0"] = n0;
    params["coulomb_log"] = coulomb_log;
    params["dt"] = dt;
    params["num_steps"] = num_steps;
    params["diag_interval"] = diag_interval;
    params["charge"] = charge;
    params["mass"] = mass;
    params["npc"] = npc;
    params["N"] = num_particles;
    params["n_density"] = n_density;
    params["T_par0"] = s0.T_par;
    params["T_perp0"] = s0.T_perp;
    std::ofstream pout(test_name + "_params.json");
    pout << params.dump(2);
    pout.close();

    Stats sf = compute_stats(sp);
    std::cout << "  Final:   T_par=" << sf.T_par << " T_perp=" << sf.T_perp << " ratio=" << sf.T_perp / sf.T_par
              << "\n";
    std::cout << "  Energy conservation: dE/E = " << (sf.Etotal - s0.Etotal) / s0.Etotal << "\n";
    std::cout << "  Output: " << fname << ", " << test_name << "_params.json\n";
}

// ============================================================================
// Test 2: two-species temperature equilibration (Te != Ti -> common T)
// ============================================================================
// To keep the distributions Maxwellian (the theory assumes it), we also run the
// intra-species collisions (e-e and i-i) in addition to the inter-species (e-i).
void run_test_equilibration(const std::string& test_name, double charge1, double mass1, double T1_0, int npc1, int np1,
                            double charge2, double mass2, double T2_0, int npc2, int np2, int num_steps, double dt,
                            double n0, double coulomb_log, int diag_interval) {
    std::cout << "\n=== " << test_name << " (equilibration) ===\n";

    auto domain_config = make_domain_config();
    Domain domain = domain_from_config(domain_config);

    auto sp1_config = make_species_config("Species1", charge1, mass1, npc1);
    auto sp2_config = make_species_config("Species2", charge2, mass2, npc2);
    ParticlesArray sp1(sp1_config, domain);
    ParticlesArray sp2(sp2_config, domain);

    const double sigma1 = std::sqrt(T1_0 / mass1);   // isotropic init
    const double sigma2 = std::sqrt(T2_0 / mass2);
    fill_maxwellian(sp1, np1, sigma1, sigma1, sigma1, 11111);
    fill_maxwellian(sp2, np2, sigma2, sigma2, sigma2, 22222);

    CoulombCollisionOperator ee("Species1", "Species1", n0, coulomb_log);
    CoulombCollisionOperator ii("Species2", "Species2", n0, coulomb_log);
    CoulombCollisionOperator ei("Species1", "Species2", n0, coulomb_log);

    Stats s0_1 = compute_stats(sp1);
    Stats s0_2 = compute_stats(sp2);
    const double E0 = s0_1.Etotal + s0_2.Etotal;

    std::cout << "  N1=" << np1 << " N2=" << np2 << " dt=" << dt << " steps=" << num_steps << "\n";
    std::cout << "  active collision channels (all 3 each step): e-e + i-i + e-i "
                 "(intra keep both species Maxwellian; inter exchanges energy)\n";
    std::cout << "  Initial: T1=" << s0_1.T_total << " T2=" << s0_2.T_total << " E=" << E0 << "\n";

    std::string fname = test_name + ".csv";
    std::ofstream out(fname);
    if (!out)
        throw std::runtime_error("Cannot write relaxation data: " + fname);
    out << std::setprecision(17);
    out << "step,time,T1,T2,Etotal,Etotal_rel_err\n";
    out << "0,0," << s0_1.T_total << "," << s0_2.T_total << "," << E0 << ",0\n";

    for (int step = 1; step <= num_steps; step++) {
        ee.collide_same_type(sp1, dt);
        ii.collide_same_type(sp2, dt);
        ei.collide_diff_type(sp1, sp2, dt);

        if (step % diag_interval == 0) {
            Stats s1 = compute_stats(sp1);
            Stats s2 = compute_stats(sp2);
            double E = s1.Etotal + s2.Etotal;
            double rel_err = (E - E0) / E0;
            out << step << "," << step * dt << "," << s1.T_total << "," << s2.T_total << "," << E << "," << rel_err
                << "\n";
        }
    }
    out.close();

    nlohmann::json params;
    params["test_type"] = "equilibration";
    params["n0"] = n0;
    params["coulomb_log"] = coulomb_log;
    params["dt"] = dt;
    params["num_steps"] = num_steps;
    params["diag_interval"] = diag_interval;
    params["charge1"] = charge1;
    params["mass1"] = mass1;
    params["npc1"] = npc1;
    params["N1"] = np1;
    params["n1_density"] = (double) np1 / npc1;
    params["T1_0"] = s0_1.T_total;
    params["charge2"] = charge2;
    params["mass2"] = mass2;
    params["npc2"] = npc2;
    params["N2"] = np2;
    params["n2_density"] = (double) np2 / npc2;
    params["T2_0"] = s0_2.T_total;
    std::ofstream pout(test_name + "_params.json");
    pout << params.dump(2);
    pout.close();

    Stats sf1 = compute_stats(sp1);
    Stats sf2 = compute_stats(sp2);
    double Ef = sf1.Etotal + sf2.Etotal;
    std::cout << "  Final:   T1=" << sf1.T_total << " T2=" << sf2.T_total << " E=" << Ef << "\n";
    std::cout << "  Energy conservation: dE/E = " << (Ef - E0) / E0 << "\n";
    std::cout << "  Output: " << fname << ", " << test_name << "_params.json\n";
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--quick") {
        // Bounded validation: resolve temperature relaxation with a fixed seed.
        // The checker supplies the independent NRL/thermal-equilibration curves.
        run_test_isotropization("quick_isotropization", -1, 1, 4e-4, 2e-4, 4000, 6, 1200, 1, 1e13, 15, 20);
        run_test_equilibration("quick_equilibration", -1, 1, 4e-4, 1, 4000, 1, 10, 4e-5, 1, 4000, 3000, 0.5, 1e13, 15,
                               50);
        return 0;
    }
    if (argc != 1) {
        std::cerr << "Usage: coulomb [--quick]\n";
        return 1;
    }
    const double n0 = 1.e13;
    const double lnL = 15.0;

    // ------------------------------------------------------------------
    // Test 1: electron temperature isotropization, T_par(0) = 2 T_perp(0).
    // Verified against the NRL / Takizuka-Abe (TA77 Eq. 20) isotropization rate:
    //   d(T_par - T_perp)/dt = -nu (T_par - T_perp).
    // npc=6 lowers the effective density so the
    // relaxation is well-resolved; N=4000 keeps statistical noise low.
    // dt chosen so <delta^2>/step ~ 0.02 << 1.
    // T in units of m_e c^2 = 511 keV: T_par0=4e-4 ~ 0.2 keV, T_perp0=2e-4 ~ 0.1 keV.
    // ------------------------------------------------------------------
    run_test_isotropization("iso_electrons", -1.0, 1.0, 4.e-4, 2.e-4, 4000, 6, 1200, 2.0, n0, lnL, 10);

    // ------------------------------------------------------------------
    // Test 2: temperature equilibration between two species, compared to the
    // Glinskiy et al. Eq.17 theory (finite mass ratio M = m_i/m_e).  All three
    // collision channels run each step (e-e, i-i, e-i) so both distributions
    // stay Maxwellian.  Reduced mass ratio mass2=100 keeps it tractable.
    // Te(0)=4e-4 (~0.2 keV), Ti(0)=4e-5 (~0.02 keV); tau_ie ~ 1.5e4 code units.
    // N=4000, npc=1 -> n_density=4000.  The residual error is a systematic TA77
    // finite-<delta^2> bias (~ dt), not statistical noise, so we run the same
    // case at two timesteps (dt=3 coarse, dt=1 refined) to show convergence to
    // theory.  Both reach t_max=15000 (~1 tau_ie); diag cadence = 150 time units.
    // ------------------------------------------------------------------
    run_test_equilibration("eq_two_species", -1.0, 1.0, 4.e-4, 1, 4000, 1.0, 100.0, 4.e-5, 1, 4000, 5000, 3.0, n0, lnL,
                           50);
    run_test_equilibration("eq_two_species_dt1", -1.0, 1.0, 4.e-4, 1, 4000, 1.0, 100.0, 4.e-5, 1, 4000, 15000, 1.0, n0,
                           lnL, 150);

    // ------------------------------------------------------------------
    // Test 3: isotropization with a stronger anisotropy (T_par/T_perp = 5).
    // Same npc=6 / dt=2 / N=4000 as Test 1.
    // ------------------------------------------------------------------
    run_test_isotropization("iso_electrons_r5", -1.0, 1.0, 5.e-4, 1.e-4, 4000, 6, 1500, 2.0, n0, lnL, 10);

    // ------------------------------------------------------------------
    // Test 4: temperature equilibration with the REAL mass ratio m_i/m_e = 1836.
    // tau_col is much larger than the mass=100 case, so we run only to partial
    // relaxation (as in TA77).  The Glinskij Eq.17 curve is overlaid for the slope.
    // dt=5, 11000 steps; smaller dt keeps the finite-<delta^2> bias low (<~0.03).
    // ------------------------------------------------------------------
    run_test_equilibration("eq_two_species_real", -1.0, 1.0, 4.e-4, 1, 4000, 1.0, 1836.0, 4.e-5, 1, 4000, 11000, 5.0,
                           n0, lnL, 100);

    std::cout << "\nAll Coulomb collision tests completed (CSV/JSON written to the current dir).\n";
    std::cout << "build.py generates the plots automatically; to redo them manually, run from\n";
    std::cout << "this dir:  python3 <srcBeren>/tests/coulomb/compare_theory.py   (theory vs sim)\n";
    std::cout << "           python3 <srcBeren>/tests/coulomb/plot_coulomb.py     (raw diagnostics)\n";
    return 0;
}
