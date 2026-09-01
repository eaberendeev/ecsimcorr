// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_H
#define SIMULATION_ECSIM_H

#include <memory>
#include <vector>

#include "DiagnosticOutput.h"
#include "ParticlesArray.h"
#include "World.h"
#include "algorithms_ecsim.h"
#include "containers.h"
#include "simulation.h"
#include "solverSLE.h"

// Main simulation class
class SimulationEcsim : public Simulation {
   public:
    SimulationEcsim(const nlohmann::json& system_config, const nlohmann::json& particles_config, int argc, char** argv)
        : Simulation(system_config, particles_config, argc, argv) {
    }
    ~SimulationEcsim();
    void init_operators() override;
    void init_fields() override;
    void prepare_step(const int timestep) override;
    void make_step(const int timestep) override;
    void make_stepNGP(const int timestep);
    virtual void diagnostic_energy(Diagnostics& diagnostic);
    void make_diagnostic(const int timestep) override;
    void prepare_block_matrix(ShapeType type);
    void convert_block_matrix(ShapeType type);
    void first_push();
    void predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& B, Field3d& J);
    void predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& E_ex, const Field3d& B, Field3d& J);
    void second_push();

    void collect_per_species_diagnostics(Diagnostics& diagnostic, double& kineticEnergy, double& kineticEnergyNew,
                                         double& totalLostEnergy);
    void compute_field_energy_and_conservation(Diagnostics& diagnostic, const IndexRange& irange, double dt,
                                               double kineticEnergy, double kineticEnergyNew, double totalLostEnergy,
                                               double totalInjectEnergy, double energyJe_ex, double dampingEnergy);

    Field3d fieldJp;        // predict current for EM solver
    Field3d fieldJp_full;   // predict current for EM solver Jp + Lmat(E+E_n);
    Field3d fieldJe;        // Esirkepov current for E correction};
    Field3d fieldE;
    Field3d fieldEn;
    Field3d fieldEp;
    Field3d fieldE_prev;   // E_{n-1}: для экстраполяции начального приближения солвера
    Field3d fieldB;
    Field3d fieldBn;
    Field3d fieldBInit;
    Field3d fieldBFull;
    Field3d fieldE_external;

    // энергия полей, поглощённая damping-слоями на текущем шаге
    double dampingEnergy_ = 0;

    // M3: предобуславливатель фиксированного оператора IMmat (строится один раз)
    std::unique_ptr<IPreconditioner> precond_;
    // M3: отдельный предобуславливатель для predictor (A = IMmat + L меняется
    // каждый шаг; по умолчанию — Jacobi из истинной диагонали A)
    std::unique_ptr<IPreconditioner> precond_pred_;

    Operator Mmat;
    Operator IMmat;
    std::vector<IndexMap> LmatX;

   private:
    std::unique_ptr<Diagnostics> diagnostic_ptr_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;
};

template <typename Func>
void for_each_particle_chess(const ParticlesArray& particles, Func&& func) {
    constexpr int CHESS_STEP = 3;
    const auto& data = particles.particlesData;
    const auto gridSize = data.size();

#pragma omp parallel
    {
        for (int xStep = 0; xStep < CHESS_STEP; ++xStep) {
            for (int yStep = 0; yStep < CHESS_STEP; ++yStep) {
#pragma omp for collapse(2) schedule(dynamic, 32)
                for (int ix = xStep; ix < gridSize.x(); ix += CHESS_STEP) {
                    for (int iy = yStep; iy < gridSize.y(); iy += CHESS_STEP) {
                        for (int iz = 0; iz < gridSize.z(); ++iz) {
                            for (auto& particle : data(ix, iy, iz)) {
                                func(particle.coord);
                            }
                        }
                    }
                }
#pragma omp barrier
            }
        }
    }
}

void update_LmatNGP(std::vector<IndexMap>& LmatX, const Vector3R& coord, const Domain& domain, double charge,
                    double mass, double mpw, const Field3d& fieldB, const double dt);
void update_Lmat(std::vector<IndexMap>& LmatX, const Vector3R& coord, const Domain& domain, double charge, double mass,
                 double mpw, const Field3d& fieldB, const double dt);
template <typename T>
void fill_matrixL(const ParticlesArray& particles, T& mat, const Field3d& fieldB, const Domain& domain, const double dt,
                  ShapeType type) {
    if (particles.is_neutral())
        return;

    const double mass = particles.mass();
    const double charge = particles.charge;
    const double mpw = particles.mpw();
    switch (type) {
        case ShapeType::NGP:
            for_each_particle_chess(particles, [&](const auto& coord) {
                update_LmatNGP(mat, coord, domain, charge, mass, mpw, fieldB, dt);
            });
            break;
        case ShapeType::Linear:
            for_each_particle_chess(
                particles, [&](const auto& coord) { update_Lmat(mat, coord, domain, charge, mass, mpw, fieldB, dt); });
            break;
        case ShapeType::Quadratic:
            std::cerr << "Fill Lmatrix for quadratic shape function is not "
                         "implemented\n";
            std::abort();
    }
}
#endif
