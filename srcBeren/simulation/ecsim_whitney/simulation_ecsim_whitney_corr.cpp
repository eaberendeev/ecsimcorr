#include "simulation_ecsim_whitney_corr.h"

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

void SimulationEcsimWhitneyCorr::init() {
    // Absorbing walls use the same clipped two-leg path as Esirkepov. Reflected
    // faces remain rejected by the base Whitney validation.
    SimulationEcsimWhitney::init();
}

void SimulationEcsimWhitneyCorr::init_fields() {
    SimulationEcsimBase::init_fields();
    fieldJe.resize(domain.size(), 3);
    fieldJe.setZero();
}

void SimulationEcsimWhitneyCorr::drift_half_step_and_assemble_response() {
    RECORD_TIMER;
    fieldBFull = fieldB + fieldBInit;
    for (auto& [name, ptr] : charged_species) {
        auto& sp = ptr.get();
        sp.currentOnGrid.setZero();
        const double qw = sp.charge() * sp.mpw();
        ParticleBoundaryEvents events;
        const PathPolicy policy = wall_boundaries_ ? bc_handler.path_policy(sp, domain, &events) : PathPolicy::none();
        for (int cell = 0; cell < sp.size(); ++cell) {
            for (auto& p : sp(cell)) {
                const Vector3R start = p.coord;
                p.move(0.5 * dt_);
                const auto path = current_particle_path(start, p, domain, policy);
                if (path.crossing.clipped) {
                    if (path.reflected)
                        throw std::runtime_error("Whitney-Corr supports absorbing walls only");
                    // Deposit only the absorbed first leg here. Survivors get
                    // their complete two-leg path in the second half.
                    const double fraction = 0.5 * path.crossing.t_hit;
                    events.record(p, {path.crossing.fate, path.crossing.face});
                    deposit_orbit_current(sp.currentOnGrid, qw, {{start, path.end, p.velocity, fraction}}, domain);
                }
            }
        }
        if (wall_boundaries_) {
            bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
            bc_handler.apply_to_particles(sp, species, domain, &events);
        } else {
            sync_periodic(sp.currentOnGrid, domain);
        }
        sp.update_cells(domain);
    }
    // The point predictor is evaluated at x^{n+1/2}; its response uses the
    // full step and the unchanged v^n, as in the ECSIM-Corr first half.
    assemble_point_response();
}

SimulationEcsimWhitneyCorr::SecondHalfStepResults
SimulationEcsimWhitneyCorr::advance_velocities_and_deposit_path_current() {
    RECORD_TIMER;
    fieldBFull = fieldB + fieldBInit;
    const Field3d electric = fieldEp + fieldE_external;
    SecondHalfStepResults result;
    for (auto& [name, ptr] : charged_species) {
        auto& sp = ptr.get();
        auto& step = result[name];
        const double qm = sp.charge() / sp.mass(), qw = sp.charge() * sp.mpw();
        const PathPolicy policy =
            wall_boundaries_ ? bc_handler.path_policy(sp, domain, &step.boundary_events) : PathPolicy::none();
        const int threads = omp_get_max_threads();
        std::vector<Field3d> currents;
        for (int t = 0; t < threads; ++t) currents.push_back(Field3d::Zero(domain.size(), 3));
        double power = 0;
#pragma omp parallel num_threads(threads) reduction(+ : power)
        {
            const int tid = omp_get_thread_num();
            Field3d& current = currents[tid];
#pragma omp for schedule(static)
            for (int cell = 0; cell < sp.size(); ++cell) {
                for (auto& p : sp(cell)) {
                    const Vector3R middle = p.coord;
                    const Vector3R old_velocity = p.velocity;
                    const Vector3R e = gather_electric(electric, middle, domain);
                    const Vector3R b = gather_magnetic(fieldBFull, middle, domain);

                    borisPusher::update_vEB(p, qm, e, b, dt_);
                    power += qw * 0.5 * (old_velocity + p.velocity).dot(e);

                    p.coord = middle + (0.5 * dt_) * p.velocity;
                    const auto path = current_particle_path(middle, p, domain, policy);
                    double fraction = 0.5;
                    Vector3R end = p.coord;
                    if (path.crossing.clipped) {
                        if (path.reflected)
                            throw std::runtime_error("Whitney-Corr supports absorbing walls only");
                        // Leave p.coord outside so the deferred boundary pass
                        // removes it with the saved wall fate.
                        end = path.end;
                        fraction = 0.5 * path.crossing.t_hit;
                        step.boundary_events.record(p, {path.crossing.fate, path.crossing.face});
                    }
                    // First-leg absorbed particles were removed before this
                    // push. Every survivor deposits its complete x^n -> x^{n+1}
                    // path here, split at the Boris turn.
                    deposit_orbit_current(current, qw,
                                          {{middle - (0.5 * dt_) * old_velocity, middle, old_velocity, 0.5},
                                           {middle, end, p.velocity, fraction}},
                                          domain);
                }
            }
        }
        step.predictor_power = power;
        for (const auto& current : currents) sp.currentOnGrid += current;
        if (wall_boundaries_)
            bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
        else
            sync_periodic(sp.currentOnGrid, domain);
    }
    return result;
}

void SimulationEcsimWhitneyCorr::solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E,
                                                                             const Field3d& B, Field3d& J, double dt) {
    RECORD_TIMER;
    Field3d rhs = E + dt * (maxwell.curlB * B - J) + field_system.M * E;
    impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);
    const double err = solve_linear_system<BicgstabSolver<Field3d>>(field_system.IM, rhs, En, E);
    LOG_STEP("  whitney corr solver error=" << err << "\n");
}

void SimulationEcsimWhitneyCorr::rescale_velocities_for_energy_balance(ParticlesArray& particles, double dt,
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
        throw std::runtime_error("invalid Whitney-Corr velocity factor for " + particles.name());
    LOG_STEP("  lambda " << particles.name() << "=" << lambda << "\n");
#pragma omp parallel for schedule(dynamic, particles.omp_granularity())
    for (int cell = 0; cell < particles.size(); ++cell)
        for (auto& p : particles(cell)) p.velocity = lambda * p.velocity;
}

void SimulationEcsimWhitneyCorr::correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step) {
    RECORD_TIMER;
    for (auto& [name, ptr] : charged_species) {
        auto& sp = ptr.get();
        const auto it = second_half_step.find(sp.name());
        if (it == second_half_step.end())
            throw std::runtime_error("missing Whitney-Corr second half for " + sp.name());
        rescale_velocities_for_energy_balance(sp, dt_, it->second.predictor_power);
        bc_handler.apply_to_particles(sp, species, domain, &it->second.boundary_events);
        sp.update_cells(domain);
    }
}

void SimulationEcsimWhitneyCorr::make_step(const int timestep) {
    RECORD_TIMER;

    continuity_.reset();
    const bool check_continuity = charge_continuity_requested(timestep);
    std::optional<Field3d> charge_before;
    if (check_continuity)
        charge_before = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());

    advance_neutral_particles(dt_);
    drift_half_step_and_assemble_response();

    {
        timer::commonTimer _t("predict_electric_field");
        predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
        bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);
    }

    auto second_half_step = advance_velocities_and_deposit_path_current();
    collect_current(fieldJe);

    {
        timer::commonTimer _t("solve_electric_field_from_current");
        solve_electric_field_from_deposited_current(fieldEn, fieldE, fieldB, fieldJe, dt_);
        bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    }

    correct_particle_energy_and_apply_boundaries(second_half_step);

    {
        timer::commonTimer _t("advance_endpoint_fields");
        maxwell.advance_magnetic_field(fieldBn, fieldB, fieldE, fieldEn, dt_);
        bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    }

    if (check_continuity) {
        const Field3d charge_after = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());
        continuity_ = measure_charge_continuity(*charge_before, charge_after, fieldJe, maxwell.divE, dt_);
    }

    bc_handler.flush_species(species);
}

void SimulationEcsimWhitneyCorr::diagnostic_energy(Diagnostics& diagnostic) const {
    diagnostic.charge_continuity = continuity_;
    measure_corr_energy(diagnostic, species, energy_fields(), bc_handler.active_range(domain.grid), dt_, dampingEnergy_,
                        filtered_current_power_);
}
