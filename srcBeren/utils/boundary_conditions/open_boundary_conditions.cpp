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

bool OpenBoundaryCondition::apply_to_particle(const Particle& p, ParticlesArray& particles,
                                              BoundaryEmitter& /*emitter*/, const Domain& domain) {
    if (!domain.geom.is_outside_face(face_, p.coord))
        return false;
    double e_kin = get_energy_particle(p.velocity, particles.mass(), particles.mpw());
    particles.diag.add_loss(face_, e_kin);
    return true;
}

bool OpenBoundaryConditionArray::apply_to_particle(const Particle& p, ParticlesArray& particles,
                                                   BoundaryEmitter& /*emitter*/, const Domain& domain) {
    for (int i = 0; i < createdBc_; ++i) {
        if (domain.geom.is_outside_face(faces_[i], p.coord)) {
            const double e_kin = get_energy_particle(p.velocity, particles.mass(), particles.mpw());
            particles.diag.add_loss(faces_[i], e_kin);
            return true;
        }
    }
    return false;
}

void OpenBoundaryConditionArray::apply_to_operator(Operator& mat, const Domain& domain) {
    // note: actually, it could write 2 times more bytes, when zeroes some elements
    RECORD_TIMER_PARAMS(sizeof(int) * (mat.rows() + mat.nonZeros()), timer::MeasureUnit::byte);

    const auto size = domain.size();
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

    const auto size = domain.size();
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
