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

    const bool doMomentum = getenvParsed<bool>("MOMENTUM_CORRECTION", false);

    if (!doMomentum) {
        // Только энергия (штатный режим)
        const double lambda = sqrt(1 + dt * (energyJe_corr - jp_cell) / energyK);

        LOG_STEP("  lambda " << sort.name() << "=" << lambda << "\n");

#pragma omp parallel for schedule(dynamic, 64)
        for (auto pk = 0; pk < sort.size(); ++pk) {
            for (auto& particle : sort.particlesData(pk)) {
                particle.velocity = lambda * particle.velocity;
            }
        }
        return;
    }

    // --- E+p коррекция (этап 2): v' = λ·v + u, замкнутая форма --------------
    // Неизвестные (λ, u): 4 скаляра на сорт; цели:
    //   K_t = K + dt·(⟨E_corr, J_e⟩ − W_push)          (энергия, как раньше)
    //   P_t = P + dt·F,  F = Σ_cells ρ·E_corr + J_e×B_avg  (импульс Лоренца)
    // Суммы: N = Σ mpw, V1 = Σ mpw·v, V2 = Σ mpw·|v|², W2 = V2 − |V1|²/N.
    // Подстановка u = (P_t/m − λ·V1)/N даёт квадратное уравнение только на λ:
    //   ½m·[λ²·W2 + |P_t|²/(m²N)] = K_t  ⇒  λ² = (K_t − |P_t|²/(2mN)) / (½m·W2).
    const double mpw = sort.mpw();
    const double m = sort.mass();
    double N = 0, V2 = 0, V1x = 0, V1y = 0, V1z = 0;
#pragma omp parallel for reduction(+ : N, V2, V1x, V1y, V1z) schedule(dynamic, 64)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        for (const auto& particle : sort.particlesData(pk)) {
            N += mpw;
            V1x += mpw * particle.velocity.x();
            V1y += mpw * particle.velocity.y();
            V1z += mpw * particle.velocity.z();
            V2 += mpw * particle.velocity.squared();
        }
    }

    // B_avg = B_n + B_ext − ¼dt·curlE·(E_n+E_{n+1})  (среднее поле шага)
    Field3d Esum = fieldE + fieldEn;
    Field3d B_avg = fieldBFull - 0.25 * dt * (mesh.curlE * Esum);

    // Плотность сорта на текущих позициях x_{n+1}: штатное обновление плотности
    // идёт ПОЗЖЕ (в particles3), и без пересчёта здесь ρ устарела на полный шаг.
    // impl начинает с densityOnGrid.setZero(), повторный вызов безопасен;
    // позиции correctv не меняет, поэтому последующий пересчёт даст то же ρ.
    sort.density_on_grid_update(SHAPE_CH);

    double Fx = 0, Fy = 0, Fz = 0;
    {
        const auto& rho = sort.densityOnGrid;
        const auto& J = sort.currentOnGrid;
        for (int i = irange.start.x(); i < irange.end.x(); ++i) {
            for (int j = irange.start.y(); j < irange.end.y(); ++j) {
                for (int k = irange.start.z(); k < irange.end.z(); ++k) {
                    const double Jx = J(i, j, k, 0), Jy = J(i, j, k, 1), Jz = J(i, j, k, 2);
                    const double Bx = B_avg(i, j, k, 0), By = B_avg(i, j, k, 1), Bz = B_avg(i, j, k, 2);
                    const double r = rho(i, j, k, 0);
                    Fx += r * fieldEp_corr_full(i, j, k, 0) + Jy * Bz - Jz * By;
                    Fy += r * fieldEp_corr_full(i, j, k, 1) + Jz * Bx - Jx * Bz;
                    Fz += r * fieldEp_corr_full(i, j, k, 2) + Jx * By - Jy * Bx;
                }
            }
        }
    }

    const double K_t = energyK + dt * (energyJe_corr - jp_cell);
    const double Ptx = m * V1x + dt * Fx, Pty = m * V1y + dt * Fy, Ptz = m * V1z + dt * Fz;
    const double W2 = V2 - (V1x * V1x + V1y * V1y + V1z * V1z) / N;
    const double num = K_t - 0.5 * (Ptx * Ptx + Pty * Pty + Ptz * Ptz) / (m * N);

    double lambda = 1.0, ux = 0.0, uy = 0.0, uz = 0.0;
    if (W2 > 1e-30 && num > 0.0) {
        lambda = sqrt(num / (0.5 * m * W2));
        ux = (Ptx / m - lambda * V1x) / N;
        uy = (Pty / m - lambda * V1y) / N;
        uz = (Ptz / m - lambda * V1z) / N;
    } else {
        // Вырождение (холодный сорт / недостижимая цель): только энергия
        lambda = sqrt(1 + dt * (energyJe_corr - jp_cell) / energyK);
    }

    LOG_STEP("  (lambda,u) " << sort.name() << " λ=" << lambda << " u=(" << ux << "," << uy << "," << uz << ")\n");

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        for (auto& particle : sort.particlesData(pk)) {
            particle.velocity = lambda * particle.velocity + Vector3R(ux, uy, uz);
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
