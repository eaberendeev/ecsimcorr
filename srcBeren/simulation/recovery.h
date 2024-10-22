// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#ifndef RECOVERY_H
#define RECOVERY_H

#include "containers.h"
#include "ParticlesArray.h"

void read_fields_from_recovery(Field3d& fieldE, Field3d& fieldB);

void write_fields_to_recovery(const Field3d& fieldE, const Field3d& fieldB,
                              const int timestep, const int recoveryInterval);

void read_particles_from_recovery(ParticlesArray& particles);
void write_particles_to_recovery(const ParticlesArray& particles,
                                 const int timestep,
                                 const int recoveryInterval);

#endif   // RECOVERY_H