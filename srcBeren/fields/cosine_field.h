#pragma once

#include <cmath>

#include "grid/domain.h"
#include "utils/containers.h"

// Adds the plane wave amplitude_c cos(k . r + phase) to component c of `field`
// at its Yee position (`type` = ELECTRIC or MAGNETIC), r from the box origin,
// k_a = 2 pi mode_a / L_a with L the box length (whole modes fit a periodic box).
inline void add_cosine_field(Field3d& field, const Domain& domain, FieldType type, const Vector3R& amplitude,
                             const Vector3R& mode, double phase) {
    const Vector3R length = domain.geom.box_max - domain.geom.box_min;
    Vector3R k;
    for (int a = 0; a < 3; ++a) k[a] = 2.0 * M_PI * mode[a] / length[a];
    const Vector3I size = field.sizes();
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int n = 0; n < size.z(); ++n)
                for (int c = 0; c < 3; ++c) {
                    const Vector3R r = domain.storage_node_position(i, j, n, type, c) - domain.geom.box_min;
                    field(i, j, n, c) += amplitude[c] * std::cos(k.dot(r) + phase);
                }
}
