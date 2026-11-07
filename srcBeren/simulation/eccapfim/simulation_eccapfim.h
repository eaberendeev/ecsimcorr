#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <map>

#include "diagnostics/Diagnostic.h"
#include "diagnostics/IDiagnosticOutput.h"
#include "simulation/ecsim/field_system.h"
#include "simulation/simulation.h"
#include "algorithms_eccapfim.h"

// x^n, v^n, E^n and B^n are integer-time states. Only an accepted nonlinear
// solution changes particles and En/Bn; iteration storage is step-local.
class SimulationEccapfim : public Simulation {
   public:
    using Simulation::Simulation;
    ~SimulationEccapfim() override;
    void init() override;
    void init_fields() override;
    void init_operators();
    // The implicit deposit uses the quadratic longitudinal shape.
    ShapeType diagnostic_density_shape() const override {
        return ShapeType::Quadratic;
    }
    void prepare_step(int timestep) override;
    void make_step(int timestep) override;
    void make_diagnostic(int timestep) override;
    void finalize() override;

    // Nodal charge density of the charged species, as deposited by the scheme.
    void collect_charge_density(Field3d& field) const;

    Field3d fieldE, fieldEn, fieldEp, fieldB, fieldBn, fieldBInit, fieldE_external, fieldJ;
    double nonlinear_residual = 0;
    int nonlinear_iterations = 0;

   private:
    using Vector = Eigen::VectorXd;
    struct ParticleIteration {
        ParticlesArray* species;
        Array3D<std::vector<Particle>> initial, trial;
        Array3D<std::vector<algorithmsECCAPFIM::OrbitBoundaryEvents>> events;
        Field3d current;
    };
    using ParticleIterations = std::map<std::string, ParticleIteration>;
    // Snapshots live only inside make_step; no iterative state is needed in recovery.
    ParticleIterations save_particles();
    void push_particles_and_deposit_current(ParticleIterations& states, const Field3d& electric,
                                            const Field3d& magnetic);
    void solve_electric_field(ParticleIterations& states);
    void accept_particles(ParticleIterations& states);
    Vector pack(const Field3d& field) const;
    void unpack(const Vector& vector, Field3d& field) const;
    IndexRange unknown_range_;
    Vector electric_free_;
    Field3d density_free_;
    EcsimFieldSystem field_system_;
    Eigen::SparseMatrix<double> maxwell_;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> maxwell_solver_;
    std::unique_ptr<Diagnostics> diagnostic_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;
    double absolute_tolerance_ = 1e-10, relative_tolerance_ = 1e-10, particle_tolerance_ = 1e-12;
    double mixing_ = 0.5;
    int max_iterations_ = 200, history_ = 6, particle_iterations_ = 100;
};
