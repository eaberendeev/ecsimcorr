#pragma once

#include "Diagnostic.h"
#include "containers.h"

// Read-only temporal levels of one completed step. Ep is the ECSIM predictor;
// prescribed fields are time independent. Species currents are prepared by the
// scheme before measuring: point current for ECSIM, path current for Corr.
struct EnergyFields {
    const Field3d& E;
    const Field3d& En;
    const Field3d& Ep;
    const Field3d& B;
    const Field3d& Bn;
    const Field3d& B_prescribed;
    const Field3d& E_external;
    const Field3d& Jp;
};

void measure_field_energy_and_conservation(Diagnostics& values, const EnergyFields& fields, const IndexRange& range,
                                           double dt, double kinetic_before, double kinetic_after, double lost_energy,
                                           double injected_energy, double external_power, double damping_energy,
                                           double filtered_current_power = 0);

void measure_ecsim_energy(Diagnostics& values, const Species& species, const EnergyFields& fields,
                          const Operator& response, const IndexRange& range, double dt, double damping_energy,
                          double filtered_current_power);
void measure_corr_energy(Diagnostics& values, const Species& species, const EnergyFields& fields,
                         const IndexRange& range, double dt, double damping_energy, double filtered_current_power);
