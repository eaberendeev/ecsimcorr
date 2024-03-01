// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H
#include <map>
#include <string>

#include "Damping.h"
#include "Diagnostic.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "Read.h"
#include "Vec.h"
#include "World.h"
#include "collision.h"

// Main simulation class
class Simulation {
   public:
    Simulation(int argc, char **argv);
    void make_all();
    void collect_charge_density(Field3d& field,
                                const std::vector<ParticlesArray>& species);

    void collect_current(Field3d& field,
                         const std::vector<ParticlesArray>& species);
};

#endif
