// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim.h"

#include "algorithms_ecsim.h"
#include "coupling/current.h"
#include "diagnostics/EnergyMeasurements.h"
#include "runtime/profiler.h"

void SimulationEcsim::drift_particles_and_assemble_response() {
    RECORD_TIMER;

    const double dt = dt_;

    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());

    for (auto& kv : species) {
        auto& sp = *kv.second;
        //  +++ x_{n-1/2} -> x_{n+1/2}
        sp.move(dt);

        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);

        // +++ get J(x_{n+1/2},v_n)_predict
        algorithmsECSIM::predict_current<SHAPE>(sp, fieldBFull, fieldJp, dt);
    }

    // Discrete filtered-current power <E0,Jp_raw>, measured before the mask.
    filtered_current_power_ = dot_product_sum(fieldE0, fieldJp, bc_handler.active_range(domain.grid));
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    assembleLmat2(dt);

    bc_handler.apply_to_operator(field_system.L, domain);
}

void SimulationEcsim::advance_particle_velocities() {
    RECORD_TIMER;

    const double dt = dt_;

    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());
    Field3d fieldE_full = fieldEp + fieldE_external;
    for (auto& kv : species) {
        auto& sp = *kv.second;
        // +++ get v'_{n+1} from v_{n} and E'_{n+1/2}
        algorithmsECSIM::predict_velocity<SHAPE>(sp, fieldE_full, fieldBFull, dt);
    }
}

void SimulationEcsim::make_step([[maybe_unused]] const int timestep) {
    RECORD_TIMER;

    drift_particles_and_assemble_response();

    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);

    advance_particle_velocities();

    // calculate fieldB
    blas::sum(2.0, fieldEp.data(), -1.0, fieldE.data(), fieldEn.data());
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);

    maxwell.advance_magnetic_field(fieldBn, fieldB, fieldE, fieldEn, dt_);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);

    // Secondaries did not contribute to this step's current/response.
    // Make them visible to endpoint diagnostics and the next step.
    bc_handler.flush_species(species);
}

void SimulationEcsim::prepare_diagnostic_currents() {
    for (auto& kv : species) {
        auto& sp = *kv.second;
        coupling::deposit_point_current(sp, sp.currentOnGrid);
        bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
    }
}

void SimulationEcsim::diagnostic_energy(Diagnostics& diagnostic) const {
    measure_ecsim_energy(diagnostic, species, energy_fields(), field_system.L, bc_handler.active_range(domain.grid),
                         dt_, dampingEnergy_, filtered_current_power_);
}
