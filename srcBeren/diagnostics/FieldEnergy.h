#pragma once

#include "containers.h"

// E and B share the code normalization: energy = 1/2 sum_active |field|^2.
// The supplied range excludes ghost copies and respects the field boundaries.
inline double calc_energy_field(const Field3d& field, const IndexRange& range) {
    return 0.5 * dot_product_sum(field, field, range);
}
