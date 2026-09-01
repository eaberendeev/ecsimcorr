#include "Diagnostic.h"
#include "log_macros.h"
#include "simulation_ecsim_corr.h"
#include "solverSLE.h"
void SimulationEcsimCorr::correctv(ParticlesArray& sort, const double dt) {
    if (sort.is_neutral())
        return;

    const Field3d fieldEp_corr_full = 0.5 * (fieldE + fieldEn) + fieldE_external;

    const auto& currentOnGrid = sort.currentOnGrid;
    const auto& domain = sort.get_domain();

    const IndexRange irange = bc_handler.active_range(domain.grid);

    const double energyJe_corr = dot_product_sum(fieldEp_corr_full, currentOnGrid, irange);

    const double jp_cell = pred_work_[sort.name()];
    const double energyK = sort.get_kinetic_energy();
    if (energyK <= 0.0)
        return;

    const double lambda = sqrt(1 + dt * (energyJe_corr - jp_cell) / energyK);

    LOG_STEP("  lambda " << sort.name() << "=" << lambda << "\n");

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        for (auto& particle : sort.particlesData(pk)) {
            particle.velocity = lambda * particle.velocity;
        }
    }
}

void SimulationEcsimCorr::correctE(Field3d& En, const Field3d& E, const Field3d& Ep, const Field3d& B, Field3d& J,
                                   const double dt) {
    RECORD_TIMER;

    timer::commonTimer timerRhs("make rhs");
    Field3d rhs = E + dt * (mesh.curlB * B - J) + mesh.Mmat * E;
    timerRhs.finish();

    // M1: начальное приближение из predictor-решения: E_{n+1} ~ 2E' - E_n
    Field3d x0 = 2.0 * Ep - E;
    if (precond_ && precond_->isExact()) {
        // One-shot Richardson с ТОЧНЫМ предобуславливателем (LU-факторизация
        // фиксированного IMmat): En = x0 + P*(rhs - IMmat*x0). Решение с
        // точностью до округления за 1 SpMV + 1 прямой подстановку.
        static const bool solverLog = getenv("SOLVER_LOG") != nullptr;
        const auto t0 = std::chrono::steady_clock::now();
        Field3d r_def(rhs.size());
        spmv(*IMmatTP_, x0, r_def);
        r_def = rhs - r_def;
        Field3d delta(rhs.size());
        precond_->apply(r_def, delta);
        const double t_apply = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        En = x0 + delta;
        if (solverLog) {
            Field3d res(rhs.size());
            res.data() = mesh.IMmat * En.data();
            const double rel = (res - rhs).norm() / rhs.norm();
            std::cout << "SOLVER [corr] iters=1 err=" << rel << " (oneshot " << t_apply << " s)\n";
        }
    } else if (precond_) {
        // M3: сильный предобуславливатель (или M2 Jacobi) для фиксированного IMmat
        solve_linear_system_precond<Field3d>(mesh.IMmat, rhs, En, x0, *precond_, "corr", IMmatTP_.get());
    } else {
        solve_linear_system<BicgstabSolver<Field3d>>(mesh.IMmat, rhs, En, x0, "corr", &mesh.IMmatDiag);
    }
    LOG_STEP("  corr solver error=" << (mesh.Imat * En - mesh.Mmat * En - rhs).norm() << "\n");
}
