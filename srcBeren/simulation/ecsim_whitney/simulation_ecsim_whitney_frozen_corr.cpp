#include "simulation_ecsim_whitney_frozen_corr.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "diagnostics/DensityMeasurements.h"
#include "diagnostics/EnergyMeasurements.h"
#include "particles/boris_pusher.h"
#include "simulation/boundary_conditions.h"
#include "runtime/log_macros.h"
#include "fields/solverSLE.h"
#include "whitney_orbit.h"
#include "whitney_shape.h"

using namespace ecsimWhitney;

void SimulationEcsimWhitneyFrozenCorr::init() {
    const auto config = system_config.value("WhitneyOrbit", nlohmann::json::object());
    if (config.value("predictor", "ballistic") != "ballistic")
        throw std::invalid_argument("ecsim_whitney_frozen_corr requires WhitneyOrbit.predictor=ballistic");
    SimulationEcsimWhitneyFrozen::init();
}

void SimulationEcsimWhitneyFrozenCorr::drift_half_step_and_assemble_response() {
    RECORD_TIMER;
    for (auto& [name, ptr] : charged_species) {
        auto& sp = ptr.get();
        // x^n -> x^{n+1/2}; the frozen response is built around this point
        // with the full-step ECSIM orbit, matching the point predictor.
        sp.move(0.5 * dt_);
        sp.update_cells(domain);
    }
    assemble_frozen_response();
}

std::map<std::string, double> SimulationEcsimWhitneyFrozenCorr::advance_velocities_and_deposit_path_current() {
    RECORD_TIMER;
    if (!response_ready_)
        throw std::runtime_error("frozen_corr velocity push requires a newly assembled response");
    const Field3d electric = fieldEp + fieldE_external;
    std::map<std::string, double> predictor_power;
    for (auto& [name, ptr] : charged_species) {
        auto& sp = ptr.get();
        sp.currentOnGrid.setZero();
        const auto& responses = response_.at(name);
        std::vector<size_t> offsets(sp.size() + 1, 0);
        for (int cell = 0; cell < sp.size(); ++cell) offsets[cell + 1] = offsets[cell] + sp.particlesData(cell).size();
        if (offsets.back() != responses.size())
            throw std::runtime_error("frozen_corr particle ordering changed before velocity push");
        const double qm = sp.charge / sp.mass(), qw = sp.charge * sp.mpw();
        const int threads = omp_get_max_threads();
        std::vector<Field3d> currents;
        for (int t = 0; t < threads; ++t) currents.push_back(Field3d::Zero(domain.size(), 3));
        auto push_cell = [&](int cell, Field3d& current) {
            double power = 0;
            size_t index = offsets[cell];
            for (auto& p : sp.particlesData(cell)) {
                const auto& response = responses[index++];
                const Vector3R middle = p.coord;
                const Vector3R old_velocity = p.velocity;
                const Vector3R e = gather_weights(electric, response.electric);
                const Vector3R b = response.magnetic;
                borisPusher::update_vEB(p, qm, e, b, dt_);
                power += qw * 0.5 * (old_velocity + p.velocity).dot(e);

                const Vector3R end = middle + (0.5 * dt_) * p.velocity;
                deposit_orbit_current(
                    current, qw,
                    {{middle - (0.5 * dt_) * old_velocity, middle, old_velocity, 0.5}, {middle, end, p.velocity, 0.5}},
                    domain);
                p.coord = end;
            }
            return power;
        };
        double power = 0;
        if (!cell_assembly_) {
            for (int cell = 0; cell < sp.size(); ++cell) power += push_cell(cell, sp.currentOnGrid);
        } else {
#pragma omp parallel num_threads(threads) reduction(+ : power)
            {
                const int tid = omp_get_thread_num();
#pragma omp for schedule(static)
                for (int cell = 0; cell < sp.size(); ++cell) power += push_cell(cell, currents[tid]);
            }
            for (const auto& current : currents) sp.currentOnGrid += current;
        }
        predictor_power[sp.name()] = power;
        sync_periodic(sp.currentOnGrid, domain);
    }
    response_ready_ = false;
    return predictor_power;
}

void SimulationEcsimWhitneyFrozenCorr::solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E,
                                                                                   const Field3d& B, Field3d& J,
                                                                                   double dt) {
    RECORD_TIMER;
    Field3d rhs = E + dt * (mesh.curlB * B - J) + mesh.Mmat * E;
    impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);
    const double err = solve_linear_system<BicgstabSolver<Field3d>>(mesh.IMmat, rhs, En, E);
    LOG_STEP("  frozen corr solver error=" << err << "\n");
}

void SimulationEcsimWhitneyFrozenCorr::rescale_velocities_for_energy_balance(ParticlesArray& particles, double dt,
                                                                             double predictor_power) {
    RECORD_TIMER;
    if (particles.is_neutral())
        return;
    const Field3d corrected_electric = 0.5 * (fieldE + fieldEn) + fieldE_external;
    const IndexRange range = bc_handler.active_range(domain.grid);
    const double corrected_power = dot_product_sum(corrected_electric, particles.currentOnGrid, range);
    const double kinetic_energy = particles.get_kinetic_energy();
    if (kinetic_energy <= 0)
        return;
    const double lambda = std::sqrt(1 + dt * (corrected_power - predictor_power) / kinetic_energy);
    if (!std::isfinite(lambda) || lambda <= 0)
        throw std::runtime_error("invalid frozen_corr velocity factor for " + particles.name());
    LOG_STEP("  lambda " << particles.name() << "=" << lambda << "\n");
#pragma omp parallel for schedule(dynamic, particles.omp_granularity())
    for (int cell = 0; cell < particles.size(); ++cell)
        for (auto& p : particles.particlesData(cell)) p.velocity = lambda * p.velocity;
}

void SimulationEcsimWhitneyFrozenCorr::correct_particle_energy_and_apply_boundaries(
    const std::map<std::string, double>& predictor_power) {
    RECORD_TIMER;
    for (auto& [name, ptr] : charged_species) {
        auto& sp = ptr.get();
        const auto it = predictor_power.find(sp.name());
        if (it == predictor_power.end())
            throw std::runtime_error("missing frozen_corr predictor power for " + sp.name());
        rescale_velocities_for_energy_balance(sp, dt_, it->second);
        bc_handler.apply_to_particles(sp, species, domain);
        sp.update_cells(domain);
    }
}

void SimulationEcsimWhitneyFrozenCorr::make_step(const int timestep) {
    RECORD_TIMER;
    globalTimer.start("Total");

    continuity_.reset();
    const bool check_continuity = charge_continuity_requested(timestep);
    std::optional<Field3d> charge_before;
    if (check_continuity)
        charge_before = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());

    advance_neutral_particles(dt_);
    drift_half_step_and_assemble_response();

    globalTimer.start("predict_electric_field");
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);
    globalTimer.finish("predict_electric_field");

    const auto predictor_power = advance_velocities_and_deposit_path_current();
    collect_current(fieldJe);

    globalTimer.start("solve_electric_field_from_current");
    solve_electric_field_from_deposited_current(fieldEn, fieldE, fieldB, fieldJe, dt_);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    globalTimer.finish("solve_electric_field_from_current");

    correct_particle_energy_and_apply_boundaries(predictor_power);

    globalTimer.start("advance_endpoint_fields");
    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, dt_);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    globalTimer.finish("advance_endpoint_fields");

    if (check_continuity) {
        const Field3d charge_after = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());
        continuity_ = measure_charge_continuity(*charge_before, charge_after, fieldJe, mesh.divE, dt_);
    }

    bc_handler.flush_species(species);
    globalTimer.finish("Total");
}

void SimulationEcsimWhitneyFrozenCorr::diagnostic_energy(Diagnostics& diagnostic) const {
    diagnostic.charge_continuity = continuity_;
    measure_corr_energy(diagnostic, species, energy_fields(), bc_handler.active_range(domain.grid), dt_, dampingEnergy_,
                        filtered_current_power_);
}
