#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "timer.h"

void PeriodicBoundaryCondition::apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) {
    RECORD_TIMER;
    if (field_type != FieldType::CURRENT && field_type != FieldType::DENSITY)
        return;

    auto sizes = field.sizes();
    auto nd = field.nd();
    const int overlap = 2 * domain.grid.ghost_cells() + 1;

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
    const auto size = domain.grid.size();
    const int overlap = 2 * domain.grid.ghost_cells() + 1;

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
