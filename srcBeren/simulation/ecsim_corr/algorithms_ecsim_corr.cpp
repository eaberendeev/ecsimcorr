#include <cmath>
#include <stdexcept>

#include "Diagnostic.h"
#include "log_macros.h"
#include "simulation_ecsim_corr.h"
#include "solverSLE.h"

void SimulationEcsimCorr::rescale_velocities_for_energy_balance(ParticlesArray& sort, const double dt,
                                                                double predictor_power) {
    RECORD_TIMER;
    if (sort.is_neutral())
        return;

    const Field3d fieldEp_corr_full = 0.5 * (fieldE + fieldEn) + fieldE_external;

    const auto& currentOnGrid = sort.currentOnGrid;
    const auto& domain = sort.get_domain();

    const IndexRange irange = bc_handler.active_range(domain.grid);

    const double corrected_field_power = dot_product_sum(fieldEp_corr_full, currentOnGrid, irange);

    const double energyK = sort.get_kinetic_energy();
    if (energyK <= 0.0)
        return;
    // K' is the kinetic energy after the Boris predictor, including particles
    // whose saved wall fate has not yet been applied. Both P values are powers
    // in code energy/time units: P_corr = <(E^n+E^{n+1})/2 + E_ext, J_Esirkepov>,
    // P_pred = sum(q*mpw*(v^n+v')/2 . E_predictor). Scaling v' by lambda gives
    // K_scaled = K' + dt*(P_corr-P_pred), hence lambda^2 = 1+dt*(P_corr-P_pred)/K'.
    // This is the discrete scheme correction; no port-work term is added.
    const double lambda = std::sqrt(1 + dt * (corrected_field_power - predictor_power) / energyK);
    if (!std::isfinite(lambda) || lambda <= 0.0)
        throw std::runtime_error("Invalid velocity correction factor for " + sort.name() +
                                 ": lambda=" + std::to_string(lambda));

    LOG_STEP("  lambda " << sort.name() << "=" << lambda << "\n");

    int64_t totalParticles = 0;
    timer::flatTimer timerWrap("OMP section wrap");
#pragma omp parallel reduction(+ : totalParticles)
    {
        timer::flatTimer timerOMP("OMP section", sort.size());
#pragma omp for schedule(dynamic, sort.omp_granularity())
        for (auto pk = 0; pk < sort.size(); ++pk) {
            std::vector<Particle>& currVec = sort.particlesData(pk);
            totalParticles += currVec.size();
            for (auto& particle : currVec) {
                particle.velocity = lambda * particle.velocity;
            }
        }

        timerOMP.m = sizeof(double) * 3 * totalParticles;
        timerOMP.unit = timer::MeasureUnit::byte;
    }
    timerWrap.m = sizeof(double) * 3 * totalParticles;
    timerWrap.unit = timer::MeasureUnit::byte;
    timerWrap.finish();
}

void SimulationEcsimCorr::solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B,
                                                                      Field3d& J, const double dt) {
    RECORD_TIMER;

    timer::commonTimer timerRhs("make rhs");
    Field3d rhs = E + dt * (maxwell.curlB * B - J) + field_system.M * E;
    // Единая Дирихле-маска er0 (mode "solv") в правой части
    impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);
    timerRhs.finish();

    const double err = solve_linear_system<BicgstabSolver<Field3d>>(field_system.IM, rhs, En, E);
    LOG_STEP("  corr solver error=" << err << "\n");
}
