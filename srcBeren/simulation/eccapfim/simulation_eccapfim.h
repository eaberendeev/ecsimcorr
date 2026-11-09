#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <map>
#include <memory>

#include <amgcl/adapter/eigen.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>

AMGCL_USE_EIGEN_VECTORS_WITH_BUILTIN_BACKEND();

#include "diagnostics/Diagnostic.h"
#include "diagnostics/IDiagnosticOutput.h"
#include "simulation/ecsim/field_system.h"
#include "simulation/simulation.h"
#include "algorithms_eccapfim.h"
#include "block_jacobi.h"

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

   protected:
    // Nonlinear source hooks. The default implementation is the kinetic
    // ECCAPFIM backend; drift-kinetic/hybrid schemes override them while
    // reusing the same Yee Maxwell operator and Anderson iteration.
    virtual void begin_nonlinear_step();
    virtual void evaluate_current(const Field3d& electric, const Field3d& magnetic);
    virtual void accept_nonlinear_step();
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
    void solve_electric_field();
    void accept_particles(ParticleIterations& states);
    ParticleIterations kinetic_states_;
    Vector pack(const Field3d& field) const;
    void unpack(const Vector& vector, Field3d& field) const;
    IndexRange unknown_range_;
    Vector electric_free_;
    Field3d density_free_;
    EcsimFieldSystem field_system_;
    Eigen::SparseMatrix<double> maxwell_;
    Eigen::SparseMatrix<double, Eigen::RowMajor> maxwell_row_;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> maxwell_solver_;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double>>
        maxwell_ic_solver_;
    using AmgBackend = amgcl::backend::builtin<double>;
    using AmgPreconditioner =
        amgcl::amg<AmgBackend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>;
    using AmgSolver = amgcl::make_solver<AmgPreconditioner, amgcl::solver::cg<AmgBackend>>;
    std::unique_ptr<AmgSolver> amgcl_solver_;
    using BlockJacobiSolver = Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                                                       eccapfim::BlockJacobi<Eigen::SparseMatrix<double>>>;
    BlockJacobiSolver maxwell_block_solver_;
    bool use_ic_ = false;
    bool use_amgcl_ = false;
    bool use_block_jacobi_ = false;
    size_t amgcl_iterations_ = 0;
    double amgcl_error_ = 0;
    std::unique_ptr<Diagnostics> diagnostic_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;
    double absolute_tolerance_ = 1e-10, relative_tolerance_ = 1e-10, particle_tolerance_ = 1e-12;
    double mixing_ = 0.5;
    int max_iterations_ = 200, history_ = 6, particle_iterations_ = 100;
    bool failure_dump_ = false, iteration_log_ = false;
    int current_timestep_ = 0;
};
