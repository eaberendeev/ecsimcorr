#pragma once

#include "config.h"
#include "containers.h"

class ParticlesArray;

namespace coupling {
// Overwrite the scalar node field with sum(weight * shape). In code density
// units, weight is q*mpw for charged species and mpw for neutral number density.
// No boundary folding is applied here; the caller applies the density BC.
// Initial and endpoint deposits must use the same charge shape as the path current.
template <ShapeType ST = SHAPE>
void deposit_density(const ParticlesArray& species, Field3d& density);

template <ShapeType ST = SHAPE>
void deposit_density_reference(const ParticlesArray& species, Field3d& density);
}   // namespace coupling
