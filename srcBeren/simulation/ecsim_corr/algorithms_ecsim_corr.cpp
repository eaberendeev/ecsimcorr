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
    if (energyK <= 0.0) return;
    const double lambda = sqrt(1 + dt * (energyJe_corr - jp_cell) / energyK);

    LOG_STEP("  lambda " << sort.name() << "=" << lambda << "\n");

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        for (auto& particle : sort.particlesData(pk)) {
            particle.velocity = lambda * particle.velocity;
        }
    }

}

void SimulationEcsimCorr::correctE(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J, const double dt) {
    RECORD_TIMER;

    timer::commonTimer timerRhs("make rhs");
    Field3d rhs = E - dt * (J + mesh.curlB * B) + mesh.Mmat * E;
    timerRhs.finish();

    solve_linear_system<BicgstabSolver<Field3d>>(mesh.IMmat, rhs, En, E);
    LOG_STEP("  corr solver error=" << (mesh.Imat * En - mesh.Mmat * En - rhs).norm() << "\n");
}
