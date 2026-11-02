#include "field_boundaries.h"

#include <numbers>

#include "simulation/boundary_conditions.h"
#include "sparse.h"
#include "timer.h"

const std::vector<ElectricDirichletEntry>& FieldBoundaryConditions::electric_dirichlet(const Domain& domain) const {
    RECORD_TIMER;

    if (electricDirichletReady_)
        return electricDirichlet_;

    electricDirichlet_.clear();
    for (const std::unique_ptr<FieldBoundaryCondition>& cond : conditions_) {
        cond->collect_electric_dirichlet(electricDirichlet_, domain);
    }

    std::stable_sort(electricDirichlet_.begin(), electricDirichlet_.end(),
                     [](const ElectricDirichletEntry& a, const ElectricDirichletEntry& b) { return a.row < b.row; });

    for (size_t i = 1; i < electricDirichlet_.size(); ++i) {
        const auto& prev = electricDirichlet_[i - 1];
        const auto& cur = electricDirichlet_[i];
        if (prev.row == cur.row && prev.value != cur.value) {
            throw std::runtime_error("conflicting prescribed electric values on DOF row " + std::to_string(cur.row) +
                                     ": " + std::to_string(prev.value) + " and " + std::to_string(cur.value));
        }
    }

    electricDirichletReady_ = true;
    return electricDirichlet_;
}

void OpenBoundaryConditionArray::apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) {
    RECORD_TIMER;

    auto set_zero = [this](double& value, int i, int j, int k, int d, const Domain& dom, FieldType ft) {
        const Vector3R pos = dom.get_node_position(i, j, k, ft, d);
        for (int numBc = 0; numBc < createdBc_; ++numBc) {
            if (dom.geom.is_outside_face(faces_[numBc], pos, -gaps_[numBc] + epss_[numBc])) {
                value = 0.0;
                break;
            }
        }
    };

    if (field_type == FieldType::CURRENT || field_type == FieldType::ELECTRIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d) set_zero(field(i, j, k, d), i, j, k, d, domain, FieldType::ELECTRIC);
    } else if (field_type == FieldType::MAGNETIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d) set_zero(field(i, j, k, d), i, j, k, d, domain, FieldType::MAGNETIC);
    } else if (field_type == FieldType::DENSITY) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k) set_zero(field(i, j, k, 0), i, j, k, 0, domain, FieldType::DENSITY);
    } else {
        std::cerr << "OpenBoundaryCondition: unsupported field type\n";
        assert(false);
    }
}
void OpenBoundaryCondition::apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) {
    RECORD_TIMER;

    auto set_zero = [gap = gap_, eps = eps_](double& value, int i, int j, int k, int d, Face face, const Domain& dom,
                                             FieldType ft) {
        const Vector3R pos = dom.get_node_position(i, j, k, ft, d);
        if (dom.geom.is_outside_face(face, pos, -gap + eps))
            value = 0.0;
    };

    if (field_type == FieldType::CURRENT || field_type == FieldType::ELECTRIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d)
                        set_zero(field(i, j, k, d), i, j, k, d, face_, domain, FieldType::ELECTRIC);
    } else if (field_type == FieldType::MAGNETIC) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    for (int d = 0; d < 3; ++d)
                        set_zero(field(i, j, k, d), i, j, k, d, face_, domain, FieldType::MAGNETIC);
    } else if (field_type == FieldType::DENSITY) {
        const Vector3I size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
        for (int i = 0; i < size.x(); ++i)
            for (int j = 0; j < size.y(); ++j)
                for (int k = 0; k < size.z(); ++k)
                    set_zero(field(i, j, k, 0), i, j, k, 0, face_, domain, FieldType::DENSITY);
    } else {
        std::cerr << "OpenBoundaryCondition: unsupported field type\n";
        assert(false);
    }
}
void OpenBoundaryConditionArray::apply_to_operator(Operator& mat, const Domain& domain) {
    // note: actually, it could write 2 times more bytes, when zeroes some elements
    RECORD_TIMER_PARAMS(sizeof(int) * (mat.rows() + mat.nonZeros()), timer::MeasureUnit::byte);

    const Vector3I size = domain.size();
    mat.makeCompressed();
    // Получаем указатели на внутренние данные CSR
    double* values = mat.valuePtr();               // Массив значений
    const int* outerIndex = mat.outerIndexPtr();   // Массив индексов начала строк
    const int* innerIndex = mat.innerIndexPtr();   // Массив индексов столбцов

#pragma omp parallel
    {
        timer::commonTimer timerOMP("omp section", 3 * (size.x() * size.y() * size.z()));
#pragma omp for schedule(dynamic, 1024)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            // Получаем диапазон ненулевых элементов для строки i
            int rowStart = outerIndex[i];
            int rowEnd = outerIndex[i + 1];
            if (rowStart == rowEnd) {
                continue;
            }

            const auto [i0, j0, k0, d0] = domain.grid.pos_vind_range(i);
            const Vector3R pos1 = domain.get_node_position(i0, j0, k0, FieldType::ELECTRIC, d0);

            bool setZeroRow = false;
            for (int bcNum = 0; bcNum < createdBc_ && !setZeroRow; ++bcNum) {
                const double eps = -gaps_[bcNum] + epss_[bcNum];
                setZeroRow = setZeroRow || domain.geom.is_outside_face(faces_[bcNum], pos1, eps);
            }

            if (setZeroRow) {
                for (int j = rowStart; j < rowEnd; j++) {
                    values[j] = 0.0;
                }
            } else {
                for (int j = rowStart; j < rowEnd; j++) {
                    int col = innerIndex[j];   // Столбец текущего элемента

                    const auto [i1, j1, k1, d1] = domain.grid.pos_vind_range(col);
                    const Vector3R pos2 = domain.get_node_position(i1, j1, k1, FieldType::ELECTRIC, d1);

                    for (int bcNum = 0; bcNum < createdBc_; ++bcNum) {
                        const double eps = -gaps_[bcNum] + epss_[bcNum];
                        const bool setZeroElem = domain.geom.is_outside_face(faces_[bcNum], pos2, eps);
                        if (setZeroElem) {
                            values[j] = 0.0;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void OpenBoundaryCondition::apply_to_operator(Operator& mat, const Domain& domain) {
    // note: actually, it could write 2 times more bytes, when zeroes some elements
    RECORD_TIMER_PARAMS(sizeof(int) * (mat.rows() + mat.nonZeros()), timer::MeasureUnit::byte);

    const Vector3I size = domain.size();
    mat.makeCompressed();
    // Получаем указатели на внутренние данные CSR
    double* values = mat.valuePtr();               // Массив значений
    const int* outerIndex = mat.outerIndexPtr();   // Массив индексов начала строк
    const int* innerIndex = mat.innerIndexPtr();   // Массив индексов столбцов

#pragma omp parallel
    {
        timer::commonTimer timerOMP("omp section", 3 * (size.x() * size.y() * size.z()));
#pragma omp for schedule(dynamic, 1024)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            // Получаем диапазон ненулевых элементов для строки i
            int rowStart = outerIndex[i];
            int rowEnd = outerIndex[i + 1];

            if (rowStart == rowEnd) {
                continue;
            }

            const auto [i0, j0, k0, d0] = domain.grid.pos_vind_range(i);

            for (int j = rowStart; j < rowEnd; j++) {
                int col = innerIndex[j];   // Столбец текущего элемента

                const auto [i1, j1, k1, d1] = domain.grid.pos_vind_range(col);
                const Vector3R pos1 = domain.get_node_position(i0, j0, k0, FieldType::ELECTRIC, d0);
                const Vector3R pos2 = domain.get_node_position(i1, j1, k1, FieldType::ELECTRIC, d1);

                const double eps = -gap_ + eps_;
                const bool setZeroRow = domain.geom.is_outside_face(face_, pos1, eps);
                const bool setZeroElem = domain.geom.is_outside_face(face_, pos2, eps);
                if (setZeroRow || setZeroElem) {
                    values[j] = 0.;
                }
            }
        }
    }
}

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
    const Vector3I size = domain.grid.size();
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

// Signed radial field in normalized units: integral Er dr = -potential_drop.
double Er0Condition::radial_field(double radius) const {
    if (profile_ == Er0Profile::Step)
        return radius >= inner_radius_ && radius <= inner_radius_ + width_ ? -potential_drop_ / width_ : 0.0;
    const double offset = radius - r0_;
    if (std::abs(offset) >= 0.5 * width_)
        return 0.0;
    return -potential_drop_ / width_ * (1.0 + std::cos(2.0 * std::numbers::pi * offset / width_));
}

void Er0Condition::collect_entries(std::vector<ElectricDirichletEntry>& out, const Domain& domain) const {
    RECORD_TIMER;

    if (entriesReady_) {
        out.insert(out.end(), entries_.begin(), entries_.end());
        return;
    }

    std::vector<ElectricDirichletEntry> built;

    const Vector3I size = domain.size();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double center_x = 0.5 * domain.num_cells().x() * dx;
    const double center_y = 0.5 * domain.num_cells().y() * dy;

    const double eps = -gap_ + eps_;
    // Кольцо задаётся только на самом граничном слое (k=ghost для ZMIN /
    // k=size_z-ghost-1 для ZMAX), остальные DOF открытого слоя получают 0.
    const int ghost = domain.grid.ghost_cells();
    const int k_face = (face_ == Face::ZMIN) ? ghost : size.z() - ghost - 1;

    for (int i = 0; i < size.x(); ++i) {
        for (int j = 0; j < size.y(); ++j) {
            for (int k = 0; k < size.z(); ++k) {
                for (int d = 0; d < 3; ++d) {
                    const Vector3R pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
                    // Только DOF открытого слоя этой грани.
                    if (!domain.geom.is_outside_face(face_, pos, eps))
                        continue;
                    // Уважаем остальные грани (цилиндр, боковые стенки).
                    if (!domain.geom.contains_ignoring_face(face_, pos, 0.0))
                        continue;

                    double value = 0.0;
                    if (d < 2 && k == k_face && width_ > 0.0) {
                        const double xx = pos.x() - center_x;
                        const double yy = pos.y() - center_y;
                        const double rr = std::hypot(xx, yy);
                        if (rr > 0.0) {
                            const double coord = (d == 0) ? xx : yy;
                            value = radial_field(rr) * coord / rr;
                        }
                    }
                    built.push_back({domain.vind(i, j, k, d), value});
                }
            }
        }
    }

    entries_ = built;
    entriesReady_ = true;
    out.insert(out.end(), entries_.begin(), entries_.end());
}

void Er0Condition::init_electric_field(Field3d& fieldE, const Domain& domain) const {
    RECORD_TIMER;

    if (face_ != Face::ZMIN && face_ != Face::ZMAX)
        return;
    if (width_ <= 0.0)
        return;

    std::vector<ElectricDirichletEntry> entries;
    collect_entries(entries, domain);

    // Пишем только ненулевые (кольцевые) значения; слой не зануляем.
    for (const auto& e : entries) {
        if (e.value == 0.0)
            continue;
        const auto [i, j, k, d] = domain.grid.pos_vind_range(e.row);
        fieldE(i, j, k, d) = e.value;
    }
}

void Er0Condition::apply_to_operator(Operator& mat, const Domain& domain) {
    RECORD_TIMER;

    if (mode_ != Er0Mode::Solve) {
        OpenBoundaryCondition::apply_to_operator(mat, domain);
        return;
    }

    mat.makeCompressed();
    double* values = mat.valuePtr();
    const int* outer = mat.outerIndexPtr();
    const int rows = 3 * domain.total_size();

    // В solv-режиме зануляем ТОЛЬКО строки открытого слоя er0-грани.
    // Столбцы не трогаем: свободные уравнения должны сохранять отклик тока
    // на предписанное граничное поле.
#pragma omp parallel for schedule(dynamic, 1024)
    for (int row = 0; row < rows; ++row) {
        const auto [i, j, k, d] = domain.grid.pos_vind_range(row);
        const Vector3R pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
        if (domain.geom.is_outside_face(face_, pos, -gap_ + eps_)) {
            for (int p = outer[row]; p < outer[row + 1]; ++p) values[p] = 0.0;
        }
    }
}
