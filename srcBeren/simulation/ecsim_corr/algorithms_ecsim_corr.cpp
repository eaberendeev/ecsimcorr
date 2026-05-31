#include "Diagnostic.h"
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
    const double lambda = sqrt(1 + dt * (energyJe_corr - jp_cell) / energyK);

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        for (auto& particle : sort.particlesData(pk)) {
            particle.velocity = lambda * particle.velocity;
        }
    }

    std::cout << "lambda " << lambda << " " << lambda * lambda << "\n";
}

void SimulationEcsimCorr::correctE(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J, const double dt) {
    Field3d rhs = E - dt * J + dt * mesh.curlB * B + mesh.Mmat * E;

    double time11 = omp_get_wtime();

    solve_linear_system<BicgstabSolver<Field3d>>(mesh.IMmat, rhs, En, E);
    double time2 = omp_get_wtime();

    std::cout << "Correction fieldE solver error = " << (mesh.Imat * En - mesh.Mmat * En - rhs).norm() << "\n";
    std::cout << "Correction fieldE Mysolver time = " << (time2 - time11) << "\n";
}
