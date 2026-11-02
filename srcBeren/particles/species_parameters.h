#pragma once

#include <string>

// Fixed for the lifetime of a species. Charge/mass and density use the code
// normalization; macro_weight = density / particles_per_cell (no cell volume).
// Restart restores particles, while these parameters come from the input.
struct SpeciesParameters {
    std::string name;
    double charge;
    double mass;
    double density;
    int particles_per_cell;
    double macro_weight;
};
