#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <map>
#include <vector>
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
#include "yee_stencil_matrix.h"

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

    Field3d fieldE, fieldEn, fieldEp, fieldB, fieldBn, fieldBHalf, fieldBInit, fieldE_external, fieldJ;
    double nonlinear_residual = 0;
    int nonlinear_iterations = 0;

   protected:
    // Nonlinear source hooks. The default implementation is the kinetic
    // ECCAPFIM backend; drift-kinetic/hybrid schemes override them while
    // reusing the same Yee Maxwell operator and Anderson iteration.
    virtual void begin_nonlinear_step();
    virtual void evaluate_current(const Field3d& electric, const Field3d& magnetic);
    virtual void accept_nonlinear_step();
    // Frozen-response mass matrix (dt/2) dJ/dE of Angus et al. (2023), with
    // particles at their latest trial positions and `magnetic` the full field
    // of the current field iterate.
    virtual void assemble_mass_matrix(const Field3d& magnetic);
    virtual double extra_kinetic_energy() const { return 0.0; }
    virtual double extra_previous_kinetic_energy() const { return 0.0; }
    // Work dt E_external . J of sources not stored in `species`.
    virtual double extra_external_work() const { return 0.0; }
    virtual void add_extra_diagnostics(Diagnostics&) {}
    int reduced_dof(int i, int j, int k, int component) const;
    // Adds W T W^T for a particle at `position`: W_c are the point weights of
    // E_c on the Yee grid (S1 along c, S2 across c), as in the orbit kernels.
    // Thread-safe: values are added atomically.
    void add_particle_mass(const Vector3R& position, const Eigen::Matrix3d& response);
    // (dt/2) M_p alone, for checks.
    Eigen::SparseMatrix<double> mass_only_matrix() const;
    void rebuild_mass_matrix(const Field3d& magnetic);
    // Point at which a particle enters the mass matrix for the orbit start -> end.
    Vector3R mass_matrix_point(const Vector3R& start, const Vector3R& end) const;
    void ensure_mass_matrix(const Field3d& magnetic);
    using Vector = Eigen::VectorXd;
    struct ParticleIteration {
        ParticlesArray* species;
        Array3D<std::vector<Particle>> initial, trial;
        Array3D<std::vector<algorithmsECCAPFIM::OrbitBoundaryEvents>> events;
        Field3d current;
        bool warm = false;   // `trial` holds orbits of a previous field iterate
    };
    using ParticleIterations = std::map<std::string, ParticleIteration>;
    // Snapshots live only inside make_step; no iterative state is needed in recovery.
    ParticleIterations save_particles();
    void push_particles_and_deposit_current(ParticleIterations& states, const Field3d& electric,
                                            const Field3d& magnetic);
    void solve_electric_field();
    void solve_fields_coupled();
    Vector pack6(const Field3d& electric, const Field3d& magnetic) const;
    void unpack6(const Vector& vector, Field3d& electric, Field3d& magnetic) const;
    Vector coupled_residual(const Vector& state);
    void solve_fields_newton_krylov(Vector& state, Vector& residual);
    void solve_fields_ngmres(Vector& state, Vector& residual);
    void solve_fields_anderson(Vector& state, Vector& residual);
    void check_mass_matrix(const Vector& state, const Vector& residual);
    Vector jacobian_vector(const Vector& state, const Vector& direction, const Vector& base_residual);
    Vector apply_right_preconditioner(const Vector& vector);
    Vector gmres(const Vector& state, const Vector& rhs, int& iterations, bool& ok);
    void accept_particles(ParticleIterations& states);
    ParticleIterations kinetic_states_;
    Vector pack(const Field3d& field) const;
    void unpack(const Vector& vector, Field3d& field) const;
    IndexRange unknown_range_;
    Vector electric_free_;
    Field3d density_free_;
    EcsimFieldSystem field_system_;
    Eigen::SparseMatrix<double> maxwell_;
    YeeStencilMatrix stencil_operator_;   // maxwell_ + (dt/2) M_p
    Eigen::SparseMatrix<double> mass_preconditioner_;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> mass_solver_;
    Vector solve_mass_system(const Vector& rhs);
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
    bool use_mass_matrix_ = false;
    bool mass_matrix_ready_ = false;
    bool use_ngmres_ = false;
    bool coupled_fields_ = false;
    // Coupled E+B nonlinear solver; the first two use the mass matrix as
    // nonlinear preconditioner.
    enum class CoupledSolver { Anderson, Ngmres, NewtonKrylov };
    CoupledSolver coupled_solver_ = CoupledSolver::NewtonKrylov;
    bool mass_matrix_every_iteration_ = false;
    bool mass_matrix_check_ = false;
    bool mass_matrix_at_end_ = false;   // xpic: latest end position; default: orbit midpoint
    bool mass_matrix_lu_ = false;   // sparse LU instead of GMRES
    std::string mass_preconditioner_name_ = "auto";
    YeeStencilMatrix::Preconditioner mass_preconditioner_kind_ = YeeStencilMatrix::Preconditioner::BlockJacobi;
    long mass_solve_iterations_ = 0;
    double mass_matrix_tolerance_ = 1e-6;
    int coupled_history_ = 30;
    int residual_evaluations_ = 0;           // particle pushes in the last step
    double coupled_atol_ = 1e-10, coupled_rtol_ = 1e-10;
    int coupled_max_outer_ = 20, coupled_gmres_restart_ = 20, coupled_gmres_max_ = 60;
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
