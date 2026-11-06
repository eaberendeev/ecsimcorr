#include <cstdlib>

#include "runtime/profiler.h"
#include "simulation/boundary_conditions.h"
#include "utils/sparse.h"

namespace {
int periodic_axis(Face face) {
    if (face == Face::XMIN || face == Face::XMAX)
        return 0;
    if (face == Face::YMIN || face == Face::YMAX)
        return 1;
    if (face == Face::ZMIN || face == Face::ZMAX)
        return 2;
    assert(false && "Periodic boundary requires a Cartesian face");
    std::abort();
}

int periodic_owner(int storage_index, int period, int ghost) {
    const int logical = storage_index - ghost;
    return ghost + (logical % period + period) % period;
}

// When the two halo bands overlap, a logical node can have more than two
// stored aliases. Sum each original deposit once before copying to aliases.
void fold_short_periodic_axis(Field3d& field, const Grid& grid, int axis) {
    const int period = grid.num_cells(axis), ghost = grid.ghost_cells();
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
                    sums[periodic_owner(k, period, ghost) - ghost] += field(node.x(), node.y(), node.z(), d);
                }
                for (int k = 0; k < sizes[axis]; ++k) {
                    node[axis] = k;
                    field(node.x(), node.y(), node.z(), d) = sums[periodic_owner(k, period, ghost) - ghost];
                }
            }
}

// Fold the deposited response onto logical owners, then broadcast each row
// to every stored alias. On a periodic E this gives the sum of the original
// response rows, including domains with several aliases in one halo band.
void fold_short_periodic_operator(Operator& mat, const Domain& domain, int axis) {
    const Grid& grid = domain.grid;
    const int period = grid.num_cells(axis), ghost = grid.ghost_cells();
    std::vector<Trip> trips;
    for (int row = 0; row < mat.outerSize(); ++row) {
        const auto [i, j, k, d] = grid.pos_vind_range(row);
        Vector3I node(i, j, k);
        const int first_alias = node[axis] % period;
        for (Operator::InnerIterator it(mat, row); it; ++it) {
            const auto [ci, cj, ck, cd] = grid.pos_vind_range(it.col());
            Vector3I column(ci, cj, ck);
            column[axis] = periodic_owner(column[axis], period, ghost);
            const int col = domain.vind(column.x(), column.y(), column.z(), cd);
            for (int alias = first_alias; alias < grid.size(axis); alias += period) {
                node[axis] = alias;
                trips.emplace_back(domain.vind(node.x(), node.y(), node.z(), d), col, it.value());
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
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
    const int overlap = 2 * domain.grid.ghost_cells() + 1;
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

void PeriodicBoundaryCondition::apply_to_operator(Operator& mat, const Domain& domain) {
    RECORD_TIMER;
    const Vector3I size = domain.grid.size();
    const int overlap = 2 * domain.grid.ghost_cells() + 1;
    const int axis = periodic_axis(face_);
    if (domain.num_cells(axis) < overlap) {
        fold_short_periodic_operator(mat, domain, axis);
        return;
    }

    Operator boundaryMatrix(mat.rows(), mat.cols());
    std::vector<Trip> boundaryTrips;
    const int last_indx = size.x() - overlap;
    const int last_indy = size.y() - overlap;
    const int last_indz = size.z() - overlap;

    if (face_ == Face::XMIN || face_ == Face::XMAX) {
        boundaryTrips.reserve(size.x() * size.y() * size.z());
#pragma omp parallel
        {
            timer::flatTimer timeOMP("OMP section X");
            std::vector<Trip> localTrips;
            localTrips.reserve(size.x() * size.y() * size.z() / omp_get_num_threads());
#pragma omp for schedule(dynamic, 32)
            for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
                auto ix = domain.pos_vind(i, 0);
                auto iy = domain.pos_vind(i, 1);
                auto iz = domain.pos_vind(i, 2);
                auto id = domain.pos_vind(i, 3);
                if (ix < overlap) {
                    auto indBound = domain.vind(last_indx + ix, iy, iz, id);

                    for (Operator::InnerIterator it(mat, i); it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = domain.pos_vind(ind2, 0);
                        auto iy1 = domain.pos_vind(ind2, 1);
                        auto iz1 = domain.pos_vind(ind2, 2);
                        auto id1 = domain.pos_vind(ind2, 3);

                        if (ix1 < overlap) {
                            auto indBound2 = domain.vind(last_indx + ix1, iy1, iz1, id1);
                            localTrips.emplace_back(indBound, indBound2, it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
                if (ix > last_indx - 1) {
                    auto indBound = domain.vind(ix - last_indx, iy, iz, id);
                    for (Operator::InnerIterator it(mat, i); it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = domain.pos_vind(ind2, 0);
                        auto iy1 = domain.pos_vind(ind2, 1);
                        auto iz1 = domain.pos_vind(ind2, 2);
                        auto id1 = domain.pos_vind(ind2, 3);
                        if (ix1 > last_indx - 1) {
                            auto indBound2 = domain.vind(ix1 - last_indx, iy1, iz1, id1);
                            localTrips.emplace_back(indBound, indBound2, it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
            }

#pragma omp critical
            boundaryTrips.insert(boundaryTrips.end(), localTrips.begin(), localTrips.end());
        }
    }
    if (face_ == Face::YMIN || face_ == Face::YMAX) {
#pragma omp parallel
        {
            timer::flatTimer timeOMP("OMP section Y");
            std::vector<Trip> localTrips;
            localTrips.reserve(size.x() * size.y() * size.z() / omp_get_num_threads());

#pragma omp for schedule(dynamic, 32)
            for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
                auto iy = domain.pos_vind(i, 1);
                auto ix = domain.pos_vind(i, 0);
                auto iz = domain.pos_vind(i, 2);
                auto id = domain.pos_vind(i, 3);
                if (iy < overlap) {
                    auto indBound = domain.vind(ix, last_indy + iy, iz, id);
                    for (Operator::InnerIterator it(mat, i); it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = domain.pos_vind(ind2, 0);
                        auto iy1 = domain.pos_vind(ind2, 1);
                        auto iz1 = domain.pos_vind(ind2, 2);
                        auto id1 = domain.pos_vind(ind2, 3);

                        if (iy1 < overlap) {
                            localTrips.emplace_back(indBound, domain.vind(ix1, last_indy + iy1, iz1, id1), it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
                if (iy > last_indy - 1) {
                    auto indBound = domain.vind(ix, iy - last_indy, iz, id);
                    for (Operator::InnerIterator it(mat, i); it; ++it) {
                        auto ind2 = it.col();
                        // auto value = it.value();
                        auto ix1 = domain.pos_vind(ind2, 0);
                        auto iy1 = domain.pos_vind(ind2, 1);
                        auto iz1 = domain.pos_vind(ind2, 2);
                        auto id1 = domain.pos_vind(ind2, 3);

                        if (iy1 > last_indy - 1) {
                            auto indBound2 = domain.vind(ix1, iy1 - last_indy, iz1, id1);
                            localTrips.emplace_back(indBound, indBound2, it.value());

                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
            }
#pragma omp critical
            boundaryTrips.insert(boundaryTrips.end(), localTrips.begin(), localTrips.end());
        }
    }
    if (face_ == Face::ZMIN || face_ == Face::ZMAX) {
#pragma omp parallel
        {
            timer::flatTimer timeOMP("OMP section Z");
            std::vector<Trip> localTrips;
            localTrips.reserve(size.x() * size.y() * size.z() / omp_get_num_threads());
#pragma omp for schedule(dynamic, 32)
            for (int i = 0; i < mat.outerSize(); i++) {
                auto iz = domain.pos_vind(i, 2);
                auto ix = domain.pos_vind(i, 0);
                auto iy = domain.pos_vind(i, 1);
                auto id = domain.pos_vind(i, 3);
                if (iz < overlap) {
                    auto indBound = domain.vind(ix, iy, last_indz + iz, id);
                    for (Operator::InnerIterator it(mat, i); it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = domain.pos_vind(ind2, 0);
                        auto iy1 = domain.pos_vind(ind2, 1);
                        auto iz1 = domain.pos_vind(ind2, 2);
                        auto id1 = domain.pos_vind(ind2, 3);

                        if (iz1 < overlap) {
                            auto indBound2 = domain.vind(ix1, iy1, last_indz + iz1, id1);
                            localTrips.emplace_back(indBound, indBound2, it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
                if (iz > last_indz - 1) {
                    auto indBound = domain.vind(ix, iy, iz - last_indz, id);
                    for (Operator::InnerIterator it(mat, i); it; ++it) {
                        auto ind2 = it.col();
                        // auto value = it.value();
                        auto ix1 = domain.pos_vind(ind2, 0);
                        auto iy1 = domain.pos_vind(ind2, 1);
                        auto iz1 = domain.pos_vind(ind2, 2);
                        auto id1 = domain.pos_vind(ind2, 3);

                        if (iz1 > last_indz - 1) {
                            auto indBound2 = domain.vind(ix1, iy1, iz1 - last_indz, id1);
                            localTrips.emplace_back(indBound, indBound2, it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
            }

#pragma omp critical
            boundaryTrips.insert(boundaryTrips.end(), localTrips.begin(), localTrips.end());
        }
    }

    timer::commonTimer timerSetFromTriplets("set from triplets");
    boundaryMatrix.setFromTriplets(boundaryTrips.begin(), boundaryTrips.end());
    timerSetFromTriplets.finish();
    mat = parallelSparseSum(mat, boundaryMatrix);
}
