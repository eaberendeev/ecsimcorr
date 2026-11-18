#pragma once

#include <array>
#include <vector>

#include "include/config.h"
#include "utils/containers.h"

class ParticlesArray;

namespace coupling {
// Overwrite the scalar node field with sum(weight * shape). In code density
// units, weight is q*mpw for charged species and mpw for neutral number density.
// No boundary folding is applied here; the caller applies the density BC.
// Initial and endpoint deposits must use the same charge shape as the path current.
// The four-node envelope starts at floor(coord/cell_size)-1 in logical nodes.
// All nonzero support must fit the configured grid halo; zero padding may extend
// beyond it. Species bins must agree with their domain grid and particle cells.
template <ShapeType ST = SHAPE>
void deposit_density(const ParticlesArray& species, Field3d& density);

// Storage cell layers per axis; a cell is selected when its index is marked on any axis.
using CellLayers = std::array<std::vector<char>, 3>;

// deposit_density from the selected cells only: exact on every node that no
// unselected cell reaches, partial elsewhere.
template <ShapeType ST = SHAPE>
void deposit_density(const ParticlesArray& species, Field3d& density, const CellLayers& cells);

template <ShapeType ST = SHAPE>
void deposit_density_reference(const ParticlesArray& species, Field3d& density);
}   // namespace coupling
