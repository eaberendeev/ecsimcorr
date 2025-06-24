// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_symplectic.h"
#include "util.h"
enum class Direction { X = 0, Y = 1, Z = 2 };

void move_particles() {
    double dt = parameters.get_double("Dt");
    for (auto& sp : species) {
        move_particles_impl<Direction::X>(sp, dt);
        move_particles_impl<Direction::Y>(sp, dt);
        move_particles_impl<Direction::Z>(sp, dt);
    }
}


template<Direction DIR>
std::pair<double, double> interpolate_magnetic_fields(int cell_idx, const auto& weights, auto& fieldB) {
    const auto& [w1, w2] = weights;
    
    if constexpr (DIR == Direction::X) {
        const int k_idx = w1.indices[0];  // wz
        const int j_idx = w2.indices[0];  // wy

        double Bz = (fieldB(cell_idx, j_idx, k_idx, Z) * w1.weights[0] +
                     fieldB(cell_idx, j_idx, k_idx + 1, Z) * w1.weights[1]);
        double By = (fieldB(cell_idx, j_idx, k_idx, Y) * w2.weights[0] +
                     fieldB(cell_idx, j_idx + 1, k_idx, Y) * w2.weights[1]);
        return {Bz, By};
        
    } else if constexpr (DIR == Direction::Y) {
        const int k_idx = w1.indices[0];  // wz
        const int i_idx = w2.indices[0];  // wx

        double Bz = (fieldB(i_idx, cell_idx, k_idx, Z) * w1.weights[0] +
                     fieldB(i_idx, cell_idx, k_idx + 1, Z) * w1.weights[1]);
        double Bx = (fieldB(i_idx, cell_idx, k_idx, X) * w2.weights[0] +
                     fieldB(i_idx + 1, cell_idx, k_idx, X) * w2.weights[1]);
        return {Bz, Bx};
        
    } else { // Direction::Z
        const int j_idx = w1.indices[0];  // wy
        const int i_idx = w2.indices[0];  // wx

        double By = (fieldB(i_idx, j_idx, cell_idx, Y) * w1.weights[0] +
                     fieldB(i_idx, j_idx + 1, cell_idx, Y) * w1.weights[1]);
        double Bx = (fieldB(i_idx, j_idx, cell_idx, X) * w2.weights[0] +
                     fieldB(i_idx + 1, j_idx, cell_idx, X) * w2.weights[1]);
        return {By, Bx};
    }
}

template <Direction DIR>
void update_current_density(int cell_idx, const auto& weights, double charge,
                            double path, double dt, auto& fieldJ) {
    const auto& [w1, w2] = weights;

#pragma omp simd collapse(2)
    for (int n1 = 0; n1 < 2; ++n1) {
        for (int n2 = 0; n2 < 2; ++n2) {
            if constexpr (DIR == Direction::X) {
                const int k_idx = w1.indices[n1];   // k_idx (z)
                const int j_idx = w2.indices[n2];   // j_idx (y)
                const double weight = w1.weights[n1] * w2.weights[n2];
#pragma omp atomic update
                fieldJ(cell_idx, j_idx, k_idx, X) +=
                    weight * charge * path / dt;

            } else if constexpr (DIR == Direction::Y) {
                const int k_idx = w1.indices[n1];   // k_idx (z)
                const int i_idx = w2.indices[n2];   // i_idx (x)
                const double weight = w1.weights[n1] * w2.weights[n2];
#pragma omp atomic update
                fieldJ(i_idx, cell_idx, k_idx, Y) +=
                    weight * charge * path / dt;

            } else {                                // Direction::Z
                const int j_idx = w1.indices[n1];   // j_idx (y)
                const int i_idx = w2.indices[n2];   // i_idx (x)
                const double weight = w1.weights[n1] * w2.weights[n2];
#pragma omp atomic update
                fieldJ(i_idx, j_idx, cell_idx, Z) +=
                    weight * charge * path / dt;
            }
        }
    }
}

template<Direction DIR>
void update_velocities(auto& velocity, double charge, double mass, 
                      double total_B1, double total_B2) {
    if constexpr (DIR == Direction::X) {
        velocity.y() -= charge / mass * total_B1;  // total_Bz
        velocity.z() += charge / mass * total_B2;  // total_By
    } else if constexpr (DIR == Direction::Y) {
        velocity.x() += charge / mass * total_B1;  // total_Bz
        velocity.z() -= charge / mass * total_B2;  // total_Bx
    } else { // Direction::Z
        velocity.x() -= charge / mass * total_B1;  // total_By
        velocity.y() += charge / mass * total_B2;  // total_Bx
    }
}

template <Direction DIR>
void SimulationSymplectic::move_particles_impl(ParticlesArray& species,
                                               const double dt) {
   // const int dir = static_cast<int>(DIR);
    const double3 cell_sizes = domain.cell_size();

    const double mass = species->mass();
    const double charge = species->charge;

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            // Получаем координаты и скорости в зависимости от направления
            const double start_pos = particle.coord(DIR);
            const double end_pos = start_pos + particle.velocity(DIR) * dt;

            const double cell_size = cell_sizes[];
            const int initial_cell =
                static_cast<int>((start_pos / cell_size) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_pos / cell_size) + GHOST_CELLS);
            const int cell_diff = final_cell - initial_cell;

            // Вычисляем веса для двух других направлений
            const auto weights =
                compute_weights_for_direction<DIR>(particle, cell_sizes);

            const auto process_cell = [&](int cell_idx, double path) {
                // Получаем магнитные поля для двух перпендикулярных направлений
                auto [B_field1, B_field2] =
                    interpolate_magnetic_fields<DIR>(cell_idx, weights, fieldB);

                // Обновляем ток в текущем направлении
                update_current_density<DIR>(cell_idx, weights, charge, path, dt,
                                            fieldJ);

                return std::make_pair(B_field1 * path, B_field2 * path);
            };

            double total_B1 = 0.0, total_B2 = 0.0;

            if (cell_diff == 0) {
                double path = end_pos - start_pos;
                auto [B1, B2] = process_cell(initial_cell, path);
                total_B1 = B1;
                total_B2 = B2;
            } else {
                const double boundary_pos =
                    cell_size * (initial_cell + (cell_diff > 0) - GHOST_CELLS);
                const double path1 = boundary_pos - start_pos;
                const double path2 = end_pos - boundary_pos;

                auto [B1_1, B2_1] = process_cell(initial_cell, path1);
                auto [B1_2, B2_2] = process_cell(final_cell, path2);

                total_B1 = B1_1 + B1_2;
                total_B2 = B2_1 + B2_2;
            }

            // Обновляем скорости в перпендикулярных направлениях
            update_velocities<DIR>(particle.velocity, charge, mass, total_B1,
                                   total_B2);

            // Обновляем координату
            particle.coord(DIR) =  end_pos;
        }
    }
}

// Обобщенная структура для хранения параметров направления
struct DirectionConfig {
    int b1_comp;
    int b2_comp;
    int vel1_comp;
    int vel2_comp;
};

// Конфигурация для каждого направления
const std::map<Direction, DirectionConfig> dir_config = {
    {Direction::X, {Z, Y, Y, Z}},
    {Direction::Y, {Z, X, Z, X}},
    {Direction::Z, {Y, X, X, Y}}};


template <Direction DIR>
void SimulationSymplectic::move_particles_impl(ParticlesArray& species,
                                               const double dt) {
    // const int dir = static_cast<int>(DIR);
    const double3 cell_sizes = domain.cell_size();

    const double mass = species->mass();
    const double charge = species->charge;
    const auto& cfg = dir_config.at(dir);
    const auto normalided_pos =
        [&](double pos, int dir) { return (pos / cell_sizes[dir]); }

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            // Получаем координаты и скорости в зависимости от направления
            const double start_pos = particle.coord(DIR);
            const double end_pos = start_pos + particle.velocity(DIR) * dt;

            const double cell_size = cell_sizes[];
            const int initial_cell =
                static_cast<int>((start_pos / cell_size) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_pos / cell_size) + GHOST_CELLS);
            const int cell_diff = final_cell - initial_cell;
            // Вычисление весов для перпендикулярных направлений
            const auto w1 = compute_weights(
                normalided_pos(coord(cfg.b1_comp) + GHOST_CELLS, cfg.b1_comp),
                0.0);
            const auto w2 = compute_weights(
                normalided_pos(coord(cfg.b2_comp) + GHOST_CELLS, cfg.b2_comp),
                0.0);
            const auto process_cell = [&](int cell_idx, double path) {
                // Получаем магнитные поля для двух перпендикулярных направлений
                auto [B_field1, B_field2] =
                    interpolate_magnetic_fields<DIR>(cell_idx, weights, fieldB);

                // Обновляем ток в текущем направлении
                update_current_density<DIR>(cell_idx, weights, charge, path, dt,
                                            fieldJ);

                return std::make_pair(B_field1 * path, B_field2 * path);
            };

            double total_B1 = 0.0, total_B2 = 0.0;

            if (cell_diff == 0) {
                double path = end_pos - start_pos;
                auto [B1, B2] = process_cell(initial_cell, path);
                total_B1 = B1;
                total_B2 = B2;
            } else {
                const double boundary_pos =
                    cell_size * (initial_cell + (cell_diff > 0) - GHOST_CELLS);
                const double path1 = boundary_pos - start_pos;
                const double path2 = end_pos - boundary_pos;

                auto [B1_1, B2_1] = process_cell(initial_cell, path1);
                auto [B1_2, B2_2] = process_cell(final_cell, path2);

                total_B1 = B1_1 + B1_2;
                total_B2 = B2_1 + B2_2;
            }

            // Обновляем скорости в перпендикулярных направлениях
            update_velocities<DIR>(particle.velocity, charge, mass, total_B1,
                                   total_B2);

            // Обновляем координату
            particle.coord(DIR) = end_pos;
        }
    }
}