#include "boundary_conditions.h"

#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "timer.h"

// -----------------------------------------------
// Конкретные реализации
// -----------------------------------------------

void PeriodicBoundaryCondition::apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) {
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
        boundaryMatrix.setFromTriplets(boundaryTrips.begin(), boundaryTrips.end());
        mat += boundaryMatrix;
        boundaryTrips.clear();
    }
    if (face_ == Face::YMIN || face_ == Face::YMAX) {
#pragma omp parallel
        {
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

        boundaryMatrix.setFromTriplets(boundaryTrips.begin(), boundaryTrips.end());
        mat += boundaryMatrix;
        boundaryTrips.clear();
    }
    if (face_ == Face::ZMIN || face_ == Face::ZMAX) {
#pragma omp parallel
        {
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

        boundaryMatrix.setFromTriplets(boundaryTrips.begin(), boundaryTrips.end());
        mat += boundaryMatrix;
        boundaryTrips.clear();
    }
}

void SecondEmissionCondition::apply_to_particle(const Particle& p, const ParticlesArray& particles,
                                                BoundaryEmitter& emitter, const Domain& domain) {
    // Условие: вторичная эмиссия работает только для электронов
    if (particles.name() != "Ions") {
        return;
    }

    if (!domain.geom.is_outside_face(face_, p.coord) && is_outside_other_faces(p.coord, domain))
        return;

    Particle p_new = p;
    p_new.velocity = gauss_.sample(pulse_gen_);

    if (face_ == Face::ZMIN) {
        p_new.coord.z() = -p_new.coord.z();

        p_new.velocity.z() = std::abs(p_new.velocity.z());
    } else if (face_ == Face::ZMAX) {
        p_new.coord.z() = 2 * domain.cell_size().z() * domain.num_cells().z() - p_new.coord.z();
        p_new.velocity.z() = std::abs(p_new.velocity.z());
    } else {
        return;
    }

    emitter.emit_to_species("Electrons", p_new);
}

void BoundaryConditionHandler::apply_to_particles(
    ParticlesArray& particles,
    [[maybe_unused]] std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species,
    const Domain& domain) {
    double time0 = omp_get_wtime();
    auto& data = particles.particlesData;
    int nx = data.size().x(), ny = data.size().y(), nz = data.size().z();

    // имя текущего сорта
    const std::string& species_name = particles.name();
    std::cout << "Applying boundary conditions to " << species_name << "\n";
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                auto& cell = data(ix, iy, iz);

                auto it = std::remove_if(cell.begin(), cell.end(), [&](Particle& p) {
                    if (domain.contains(p.coord))
                        return false;

                    if (periodic_[0] || periodic_[1] || periodic_[2]) {
                        Particle p_new = p;
                        p_new.coord = wrap_periodic(p.coord, domain);

                        if (domain.contains(p_new.coord))
                            emitter.emit_current_species(p_new);
                    }

                    for (const auto& cond : conditions_) {
                        cond->apply_to_particle(p, particles, emitter, domain);
                    }
                    return true;
                });

                cell.erase(it, cell.end());
            }
        }
    }
    int count = 0;
    // добавить новые частицы текущего сорта
    for (const auto& p : emitter.current_species_particles()) {
        particles.add_particle(p);
        count++;
    }
    emitter.clear_current_species_buffer();
    std::cout << "Added " << count << " " << species_name << " particles from emitter\n";
    std::cout << "Time to apply boundary conditions: " << omp_get_wtime() - time0 << " seconds\n";
}

void BoundaryConditionHandler::flush_species(
    std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) {
    // добавить частицы в другие сорта по имени
    const auto& other = emitter.other_species_particles();

    for (const auto& kv : other) {
        int count = 0;
        const std::string& name = kv.first;
        if (all_species.count(name) == 0) {
            throw std::runtime_error("Species " + name + " not found when adding from emitter");
        }
        auto& target_species = all_species[name];
        for (const auto& p : kv.second) {
            target_species->add_particle(p);
            count++;
        }
        std::cout << "Added " << count << " " << name << " particles from emitter\n";
    }
    emitter.clear_other_species_buffers();
}

void OpenBoundaryCondition::apply_to_operator(Operator& mat, const Domain& domain) {
    RECORD_TIMER;

    const auto size = domain.size();
    mat.makeCompressed();
    // Получаем указатели на внутренние данные CSR
    double* values = mat.valuePtr();               // Массив значений
    const int* outerIndex = mat.outerIndexPtr();   // Массив индексов начала строк
    const int* innerIndex = mat.innerIndexPtr();   // Массив индексов столбцов

    auto set_values_zero = [gap = gap_, eps = eps_](double* values, Face face, const Domain& domain, int vindg,
                                                    int vindg1, int value_index) {
        int i = domain.grid.pos_vind(vindg, 0);
        int j = domain.grid.pos_vind(vindg, 1);
        int k = domain.grid.pos_vind(vindg, 2);
        int d = domain.grid.pos_vind(vindg, 3);
        int i1 = domain.grid.pos_vind(vindg1, 0);
        int j1 = domain.grid.pos_vind(vindg1, 1);
        int k1 = domain.grid.pos_vind(vindg1, 2);
        int d1 = domain.grid.pos_vind(vindg1, 3);
        auto pos1 = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
        auto pos2 = domain.get_node_position(i1, j1, k1, FieldType::ELECTRIC, d1);
        bool setZero1 = domain.geom.is_outside_face(face, pos1, -gap + eps);
        bool setZero2 = domain.geom.is_outside_face(face, pos2, -gap + eps);
        if (setZero1 || setZero2) {
            values[value_index] = 0.;
        }
    };

#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
        // Получаем диапазон ненулевых элементов для строки i
        int rowStart = outerIndex[i];
        int rowEnd = outerIndex[i + 1];

        for (int j = rowStart; j < rowEnd; j++) {
            int col = innerIndex[j];   // Столбец текущего элемента
            set_values_zero(values, face_, domain, i, col, j);
        }
    }
}

void BphiCondition::apply_to_particle(const Particle& p, const ParticlesArray& particles, BoundaryEmitter& emitter,
                                      const Domain& domain) {
    // Определяем, вышла ли частица через нужную Z-грань
    if (face_ != Face::ZMIN && face_ != Face::ZMAX)
        return;   // на всякий случай

    const double eps = 1e-12;

    // Проверяем, что точка действительно вне области по Z (с учётом возможного
    // gap_)
    if (!domain.geom.is_outside_face(face_, p.coord, eps))
        return;

    // Дополнительная проверка: точка должна быть внутри всех остальных граней
    // (т.е. не выходить также по X/Y/цилиндру)
    if (!domain.geom.contains_ignoring_face(face_, p.coord, eps))
        return;

    double vz = p.velocity.z();

    bool inside_circle = is_inside_central_circle(p.coord, domain);
    double mass = particles.mass();
    bool reflect = should_electron_reflect(inside_circle, vz, mass);

    // Ветвление по сорту частиц (электроны/ионы) как в оригинале
    if (particles.name() == "Electrons" && reflect) {
        Particle new_p = p;
        new_p.velocity.z() = -vz;
        new_p.coord = domain.geom.reflect_from_face(face_, p.coord);
        // Добавляем обратно в тот же сорт через эмиттер
        emitter.emit_current_species(new_p);
        return;
    }

    // Индексы ячейки для накопления тока
    auto idxs = domain.grid.get_field_node_index(p.coord, FieldType::CURRENT, Z);
    int i = idxs.x();
    int j = idxs.y();
    // Учитываем заряд и вес частицы
    double charge = particles.charge;
    double mpw = particles.mpw();
    double dJ = charge * mpw * vz;

    Jz_(i, j) += dJ;
    return;
}

void BphiCondition::set_Bphi(Field3d& fieldB, const Array2D<double>& Jz, int k, const Domain& domain) {
    const auto size_x = fieldB.sizes().x();
    const auto size_y = fieldB.sizes().y();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            double Bx = 0;
            double By = 0;
            for (int i1 = 0; i1 < size_x; ++i1) {
                for (int j1 = 0; j1 < size_y; ++j1) {
                    double x = i * dx;
                    double y = (j + 0.5) * dy;
                    double x1 = i1 * dx;
                    double y1 = j1 * dy;

                    double rr2 = (x - x1) * (x - x1) + (y - y1) * (y - y1);
                    Bx += -Jz(i1, j1) * (y - y1) * dx * dy / (2 * M_PI * rr2);
                }
            }
            for (int i1 = 0; i1 < size_x; ++i1) {
                for (int j1 = 0; j1 < size_y; ++j1) {
                    double x = (i + 0.5) * dx;
                    double y = j * dy;
                    double x1 = i1 * dx;
                    double y1 = j1 * dy;
                    double rr2 = (x - x1) * (x - x1) + (y - y1) * (y - y1);
                    By += Jz(i1, j1) * (x - x1) * dx * dy / (2 * M_PI * rr2);
                }
            }
            fieldB(i, j, k, Axis::X) = Bx;
            fieldB(i, j, k, Axis::Y) = By;
        }
    }
}
