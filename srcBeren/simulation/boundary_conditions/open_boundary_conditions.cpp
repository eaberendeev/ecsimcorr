#include "simulation/boundary_conditions.h"

#include "runtime/profiler.h"

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

ParticleFateResult OpenBoundaryConditionArray::classify_particle(const Particle& p, const ParticlesArray& /*particles*/,
                                                                 const Domain& domain, Face hit_face) const {
    for (int i = 0; i < createdBc_; ++i) {
        if (hit_face != Face::UNAVAILABLE ? hit_face == faces_[i] : domain.geom.is_outside_face(faces_[i], p.coord)) {
            return {ParticleFate::Removed, faces_[i]};
        }
    }
    return {};
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

ParticleFateResult OpenBoundaryCondition::classify_particle(const Particle& p, const ParticlesArray& /*particles*/,
                                                            const Domain& domain, Face hit_face) const {
    if (hit_face != Face::UNAVAILABLE ? hit_face != face_ : !domain.geom.is_outside_face(face_, p.coord))
        return {};
    return {ParticleFate::Removed, face_};
}
