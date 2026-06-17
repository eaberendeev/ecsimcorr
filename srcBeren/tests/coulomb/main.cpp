#include <omp.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "ParticlesArray.h"
#include "World.h"
#include "collision.h"
#include "collision_manager.h"

struct Stats {
    double Tx, Ty, Tz;
    double Etotal;
};

Stats compute_stats(ParticlesArray &sp) {
    double sum_vx2 = 0, sum_vy2 = 0, sum_vz2 = 0;
    int count = 0;
    for (int pk = 0; pk < sp.size(); pk++) {
        for (const auto &p : sp.particlesData(pk)) {
            sum_vx2 += p.velocity.x() * p.velocity.x();
            sum_vy2 += p.velocity.y() * p.velocity.y();
            sum_vz2 += p.velocity.z() * p.velocity.z();
            count++;
        }
    }
    if (count == 0) return {0, 0, 0, 0};
    Stats s;
    s.Tx = sum_vx2 / count;
    s.Ty = sum_vy2 / count;
    s.Tz = sum_vz2 / count;
    s.Etotal = 0.5 * sp.mass() * (sum_vx2 + sum_vy2 + sum_vz2);
    return s;
}

Stats compute_stats_two(ParticlesArray &sp1, ParticlesArray &sp2) {
    double E1 = 0, E2 = 0;
    double T1x = 0, T1y = 0, T1z = 0;
    double T2x = 0, T2y = 0, T2z = 0;
    int c1 = 0, c2 = 0;
    for (int pk = 0; pk < sp1.size(); pk++) {
        for (const auto &p : sp1.particlesData(pk)) {
            T1x += p.velocity.x() * p.velocity.x();
            T1y += p.velocity.y() * p.velocity.y();
            T1z += p.velocity.z() * p.velocity.z();
            c1++;
        }
    }
    for (int pk = 0; pk < sp2.size(); pk++) {
        for (const auto &p : sp2.particlesData(pk)) {
            T2x += p.velocity.x() * p.velocity.x();
            T2y += p.velocity.y() * p.velocity.y();
            T2z += p.velocity.z() * p.velocity.z();
            c2++;
        }
    }
    Stats s;
    s.Tx = (c1 > 0) ? T1x / c1 : 0;
    s.Ty = (c2 > 0) ? T2x / c2 : 0;
    s.Tz = 0;
    E1 = 0.5 * sp1.mass() * (T1x + T1y + T1z);
    E2 = 0.5 * sp2.mass() * (T2x + T2y + T2z);
    s.Etotal = E1 + E2;
    return s;
}

void fill_maxwellian(ParticlesArray &sp, int cell, int num_particles, double sigma_x, double sigma_y, double sigma_z,
                     unsigned seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<double> dist_x(0.0, sigma_x);
    std::normal_distribution<double> dist_y(0.0, sigma_y);
    std::normal_distribution<double> dist_z(0.0, sigma_z);

    for (int i = 0; i < num_particles; i++) {
        Particle p;
        p.coord = Vector3R(0.5, 0.5, 0.5);
        p.velocity = Vector3R(dist_x(gen), dist_y(gen), dist_z(gen));
        sp.particlesData(cell).push_back(p);
    }
}

nlohmann::json make_domain_config() {
    nlohmann::json config;
    config["Dx"] = 1.0;
    config["Dy"] = 1.0;
    config["Dz"] = 1.0;
    config["NumCellsX"] = 2;
    config["NumCellsY"] = 2;
    config["NumCellsZ"] = 2;
    return config;
}

nlohmann::json make_species_config(const std::string &name, double charge, double mass, int npc) {
    nlohmann::json config;
    config["Name"] = name;
    config["Charge"] = charge;
    config["Mass"] = mass;
    config["Density"] = 1.0;
    config["NumPartPerCell"] = npc;
    return config;
}

void run_test_same_type(const std::string &test_name, double charge, double mass, double sigma_x, double sigma_yz,
                        int num_particles, int npc, int num_steps, double dt, double n0, int diag_interval) {
    std::cout << "\n=== " << test_name << " ===\n";

    auto domain_config = make_domain_config();
    Domain domain;
    domain.init_from_json(domain_config);

    auto sp_config = make_species_config("TestSpecies", charge, mass, npc);
    ParticlesArray sp(sp_config, domain);

    fill_maxwellian(sp, 0, num_particles, sigma_x, sigma_yz, sigma_yz, 12345);

    CoulombCollisionOperator collider("TestSpecies", "TestSpecies", n0, 15.0);

    Stats s0 = compute_stats(sp);
    std::cout << "  n_eff=" << (double)num_particles / npc << " dt=" << dt << " steps=" << num_steps << "\n";
    std::cout << "  Initial: Tx=" << s0.Tx << " Ty=" << s0.Ty << " Tz=" << s0.Tz << " E=" << s0.Etotal << "\n";

    std::string fname = test_name + ".csv";
    std::ofstream out(fname);
    out << "step,Tx,Ty,Tz,Etotal,Etotal_rel_err\n";
    out << "0," << s0.Tx << "," << s0.Ty << "," << s0.Tz << "," << s0.Etotal << ",0\n";

    for (int step = 1; step <= num_steps; step++) {
        collider.collide_same_type(sp, dt);

        if (step % diag_interval == 0) {
            Stats s = compute_stats(sp);
            double rel_err = (s.Etotal - s0.Etotal) / s0.Etotal;
            out << step << "," << s.Tx << "," << s.Ty << "," << s.Tz << "," << s.Etotal << "," << rel_err << "\n";
        }
    }
    out.close();

    Stats sf = compute_stats(sp);
    double rel_err = (sf.Etotal - s0.Etotal) / s0.Etotal;
    std::cout << "  Final:   Tx=" << sf.Tx << " Ty=" << sf.Ty << " Tz=" << sf.Tz << " E=" << sf.Etotal << "\n";
    std::cout << "  Energy conservation: dE/E = " << rel_err << "\n";
    double T_avg = (sf.Tx + sf.Ty + sf.Tz) / 3.0;
    double anisotropy = std::max({std::abs(sf.Tx - T_avg), std::abs(sf.Ty - T_avg), std::abs(sf.Tz - T_avg)}) / T_avg;
    std::cout << "  Anisotropy: " << anisotropy * 100 << "%\n";
    std::cout << "  Output: " << fname << "\n";
}

void run_test_diff_type(const std::string &test_name, double charge1, double mass1, double sigma1, int npc1,
                        double charge2, double mass2, double sigma2, int npc2, int np1, int np2, int num_steps,
                        double dt, double n0, int diag_interval) {
    std::cout << "\n=== " << test_name << " ===\n";

    auto domain_config = make_domain_config();
    Domain domain;
    domain.init_from_json(domain_config);

    auto sp1_config = make_species_config("Species1", charge1, mass1, npc1);
    auto sp2_config = make_species_config("Species2", charge2, mass2, npc2);
    ParticlesArray sp1(sp1_config, domain);
    ParticlesArray sp2(sp2_config, domain);

    fill_maxwellian(sp1, 0, np1, sigma1, sigma1, sigma1, 11111);
    fill_maxwellian(sp2, 0, np2, sigma2, sigma2, sigma2, 22222);

    CoulombCollisionOperator collider("Species1", "Species2", n0, 15.0);

    Stats s0_1 = compute_stats(sp1);
    Stats s0_2 = compute_stats(sp2);
    double E0 = s0_1.Etotal + s0_2.Etotal;
    double T1_0 = (s0_1.Tx + s0_1.Ty + s0_1.Tz) / 3.0;
    double T2_0 = (s0_2.Tx + s0_2.Ty + s0_2.Tz) / 3.0;

    std::cout << "  n1_eff=" << (double)np1 / npc1 << " n2_eff=" << (double)np2 / npc2 << " dt=" << dt << "\n";
    std::cout << "  Initial: T1=" << T1_0 << " T2=" << T2_0 << " E=" << E0 << "\n";

    std::string fname = test_name + ".csv";
    std::ofstream out(fname);
    out << "step,T1,T2,Etotal,Etotal_rel_err\n";
    out << "0," << T1_0 << "," << T2_0 << "," << E0 << ",0\n";

    for (int step = 1; step <= num_steps; step++) {
        collider.collide_diff_type(sp1, sp2, dt);

        if (step % diag_interval == 0) {
            Stats s1 = compute_stats(sp1);
            Stats s2 = compute_stats(sp2);
            double E = s1.Etotal + s2.Etotal;
            double T1 = (s1.Tx + s1.Ty + s1.Tz) / 3.0;
            double T2 = (s2.Tx + s2.Ty + s2.Tz) / 3.0;
            double rel_err = (E - E0) / E0;
            out << step << "," << T1 << "," << T2 << "," << E << "," << rel_err << "\n";
        }
    }
    out.close();

    Stats sf1 = compute_stats(sp1);
    Stats sf2 = compute_stats(sp2);
    double Ef = sf1.Etotal + sf2.Etotal;
    double T1_f = (sf1.Tx + sf1.Ty + sf1.Tz) / 3.0;
    double T2_f = (sf2.Tx + sf2.Ty + sf2.Tz) / 3.0;
    std::cout << "  Final:   T1=" << T1_f << " T2=" << T2_f << " E=" << Ef << "\n";
    std::cout << "  Energy conservation: dE/E = " << (Ef - E0) / E0 << "\n";
    std::cout << "  Output: " << fname << "\n";
}

void run_test_old_same_type(const std::string &test_name, double charge, double mass, double sigma_x, double sigma_yz,
                            int num_particles, int npc, int num_steps, double dt, double n0, int diag_interval) {
    std::cout << "\n=== " << test_name << " (old impl) ===\n";

    auto domain_config = make_domain_config();
    Domain domain;
    domain.init_from_json(domain_config);

    auto sp_config = make_species_config("TestSpecies", charge, mass, npc);
    auto sp_ptr = std::make_unique<ParticlesArray>(sp_config, domain);

    fill_maxwellian(*sp_ptr, 0, num_particles, sigma_x, sigma_yz, sigma_yz, 12345);

    BinaryCollider collider(n0);

    Stats s0 = compute_stats(*sp_ptr);
    std::cout << "  Initial: Tx=" << s0.Tx << " Ty=" << s0.Ty << " Tz=" << s0.Tz << " E=" << s0.Etotal << "\n";

    std::string fname = test_name + "_old.csv";
    std::ofstream out(fname);
    out << "step,Tx,Ty,Tz,Etotal,Etotal_rel_err\n";
    out << "0," << s0.Tx << "," << s0.Ty << "," << s0.Tz << "," << s0.Etotal << ",0\n";

    Species species_map;
    species_map.emplace("TestSpecies", std::move(sp_ptr));
    auto &sp = *species_map["TestSpecies"];

    for (int step = 1; step <= num_steps; step++) {
        collider.collide_same_sort_binary(species_map, dt);

        if (step % diag_interval == 0) {
            Stats s = compute_stats(sp);
            double rel_err = (s.Etotal - s0.Etotal) / s0.Etotal;
            out << step << "," << s.Tx << "," << s.Ty << "," << s.Tz << "," << s.Etotal << "," << rel_err << "\n";
        }
    }
    out.close();

    Stats sf = compute_stats(sp);
    double rel_err = (sf.Etotal - s0.Etotal) / s0.Etotal;
    std::cout << "  Final:   Tx=" << sf.Tx << " Ty=" << sf.Ty << " Tz=" << sf.Tz << " E=" << sf.Etotal << "\n";
    std::cout << "  Energy conservation: dE/E = " << rel_err << "\n";
    std::cout << "  Output: " << fname << "\n";
}

int main() {
    const double n0 = 1.e13;
    const int N = 1000;

    // Test 1: e-e thermalization (anisotropic -> isotropic)
    // npc=100 → n_eff=10, σ_x=0.01, dt=100 → variance~1.5/step → fast thermalization
    run_test_same_type("ee_thermalization", -1.0, 1.0, 0.01, 0.003, N, 100, 2000, 100.0, n0, 5);

    // Same test with old implementation for reference
    run_test_old_same_type("ee_thermalization", -1.0, 1.0, 0.01, 0.003, N, 100, 2000, 100.0, n0, 5);

    // Test 2: i-i thermalization
    // μ_ii = m_i/2 = 918 → need n_eff=1000 to compensate m² factor
    // npc=1 → n_eff=1000, σ=0.001, dt=1000
    run_test_same_type("ii_thermalization", 1.0, 1836.0, 0.001, 0.0003, N, 1, 500, 1000.0, n0, 2);

    // Test 3: e-i energy exchange
    // μ_ei ≈ 1, so collision rate similar to e-e.
    // But energy exchange rate ∝ me/mi per collision → slow transfer.
    // npc=50 → n_eff=20, σ_e=0.008, σ_i=0.0003, dt=200
    run_test_diff_type("ei_exchange", -1.0, 1.0, 0.008, 50, 1.0, 1836.0, 0.0003, 50, N, N, 5000, 200.0, n0, 10);

    std::cout << "\nAll Coulomb collision tests completed.\n";
    std::cout << "Run plot_coulomb.py to visualize results.\n";
    return 0;
}
