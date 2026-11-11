#pragma once

#include <map>
#include <string>

#include "coupling/shape.h"
#include "utils/containers.h"

// Grid moments of a species that is not stored as a ParticlesArray (the
// guiding centres of the drift-kinetic model). The definitions follow the
// kinetic diagnostics, so outputs can treat both kinds of species alike.
class SpeciesMomentSource {
   public:
    virtual ~SpeciesMomentSource() = default;

    virtual bool is_neutral() const = 0;
    // Charge density q w S(x) of the markers at storage nodes with the
    // diagnostic shape (Linear or Quadratic), overwriting `density`. The caller
    // applies the DENSITY boundary condition, as for kinetic species.
    virtual void deposit_density(Field3d& density, ShapeType shape) const = 0;
    // Current density of the species on the E positions of the Yee grid, at
    // the time level of the deposited kinetic current (the accepted step).
    virtual const Field3d& current() const = 0;
    // Raw second velocity moment m w <v_a v_b> with the linear scalar-node
    // shape, overwriting `moment`; a, b = 0, 1, 2 are the cylindrical r, phi, z
    // components (the same basis as the kinetic second moment).
    virtual void deposit_second_moment(Field3d& moment, int a, int b) const = 0;
};

using MomentSources = std::map<std::string, const SpeciesMomentSource*>;
