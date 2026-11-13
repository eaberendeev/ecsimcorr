#include "runtime/profiler.h"
#include "simulation/boundary_conditions.h"

namespace {
int periodic_axis(Face face) {
    if (face == Face::XMIN || face == Face::XMAX)
        return 0;
    if (face == Face::YMIN || face == Face::YMAX)
        return 1;
    if (face == Face::ZMIN || face == Face::ZMAX)
        return 2;
    assert(false && "Periodic boundary requires a Cartesian face");
    return 0;
}

// When the two halo bands overlap, a logical node can have more than two
// stored aliases. Sum each original deposit once before copying to aliases.
void fold_short_periodic_axis(Field3d& field, const Grid& grid, int axis) {
    const int period = grid.num_cells(axis), ghost = GHOST_CELLS;
    const Vector3I sizes = field.sizes();
    const int a = (axis + 1) % 3, b = (axis + 2) % 3;
    std::vector<double> sums(period);
    for (int i = 0; i < sizes[a]; ++i)
        for (int j = 0; j < sizes[b]; ++j)
            for (int d = 0; d < field.nd(); ++d) {
                std::fill(sums.begin(), sums.end(), 0.0);
                Vector3I node;
                node[a] = i;
                node[b] = j;
                for (int k = 0; k < sizes[axis]; ++k) {
                    node[axis] = k;
                    sums[grid.periodic_owner(k, axis) - ghost] += field(node.x(), node.y(), node.z(), d);
                }
                for (int k = 0; k < sizes[axis]; ++k) {
                    node[axis] = k;
                    field(node.x(), node.y(), node.z(), d) = sums[grid.periodic_owner(k, axis) - ghost];
                }
            }
}

}   // namespace

void PeriodicBoundaryCondition::apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) {
    RECORD_TIMER;
    if (field_type != FieldType::CURRENT && field_type != FieldType::DENSITY)
        return;

    auto sizes = field.sizes();
    auto nd = field.nd();
    // Storage has N physical periods plus the duplicated endpoint and both
    // halos: size=N+2*g+1. Opposite alias bands are separated by exactly N.
    const int overlap = 2 * GHOST_CELLS + 1;
    const int axis = periodic_axis(face_);
    if (domain.num_cells(axis) < overlap) {
        fold_short_periodic_axis(field, domain.grid, axis);
        return;
    }

    if (face_ == Face::XMIN || face_ == Face::XMAX) {
        auto i_max = sizes.x() - overlap;
        for (auto i = 0; i < overlap; ++i) {
            for (auto j = 0; j < sizes.y(); ++j) {
                for (auto k = 0; k < sizes.z(); ++k) {
                    for (auto dim = 0; dim < nd; dim++) {
                        field(i, j, k, dim) += field(i + i_max, j, k, dim);
                        field(i + i_max, j, k, dim) = field(i, j, k, dim);
                    }
                }
            }
        }
    }
    if (face_ == Face::YMIN || face_ == Face::YMAX) {
        auto j_max = sizes.y() - overlap;
        for (auto i = 0; i < sizes.x(); ++i) {
            for (auto j = 0; j < overlap; ++j) {
                for (auto k = 0; k < sizes.z(); ++k) {
                    for (auto dim = 0; dim < nd; dim++) {
                        field(i, j, k, dim) += field(i, j + j_max, k, dim);
                        field(i, j + j_max, k, dim) = field(i, j, k, dim);
                    }
                }
            }
        }
    }
    if (face_ == Face::ZMIN || face_ == Face::ZMAX) {
        auto k_max = sizes.z() - overlap;
        for (auto i = 0; i < sizes.x(); ++i) {
            for (auto j = 0; j < sizes.y(); ++j) {
                for (auto k = 0; k < overlap; ++k) {
                    for (auto dim = 0; dim < nd; dim++) {
                        field(i, j, k, dim) += field(i, j, k + k_max, dim);
                        field(i, j, k + k_max, dim) = field(i, j, k, dim);
                    }
                }
            }
        }
    }
}
