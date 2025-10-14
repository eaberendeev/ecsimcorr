// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_symplectic.h"

void SimulationSymplectic::move_x(std::unique_ptr<ParticlesArray>& species,
                                  const Field3d& fieldB, Field3d& fieldJ,
                                  const double dt) {
    if (species->is_neutral())
        return;
    double xCellSize = domain.cell_size().x();
    double yCellSize = domain.cell_size().y();
    double zCellSize = domain.cell_size().z();
    double mass = species->mass();
    double charge = species->charge;
    double mpw = species->mpw();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            const double start_x = particle.coord.x();
            const double end_x = start_x + particle.velocity.x() * dt;
            const int initial_cell =
                static_cast<int>((start_x / xCellSize) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_x / xCellSize) + GHOST_CELLS);
            const int cell_diff = final_cell - initial_cell;


            const auto wy = compute_weights(
                (particle.coord.y() / yCellSize) + GHOST_CELLS, 0.0);
            const auto wz = compute_weights(
                (particle.coord.z() / zCellSize) + GHOST_CELLS, 0.0);

            const auto process_cell = [&](int cell_idx, double path) {
                const int j_idx = wy.indices[0];
                const int k_idx = wz.indices[0];
                double Bz =
                    (fieldB(cell_idx, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(cell_idx, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By =
                    (fieldB(cell_idx, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(cell_idx, j_idx + 1, k_idx, Y) * wy.weights[1]);

#pragma omp simd collapse(2)
                for (int ny = 0; ny < 2; ++ny) {
                    for (int nz = 0; nz < 2; ++nz) {
                        const int j_idx = wy.indices[ny];
                        const int k_idx = wz.indices[nz];
                        // if (cell_idx < 0 || cell_idx >= 23 ||
                        // j_idx < 0 || j_idx >= 23 ||
                        // k_idx < 0 || k_idx >= 23) {
                        //     std::cout
                        //         << "Particle out of bounds: " << final_cell
                        //         << std::endl;
                        //     std::cout << particle << std::endl;
                        // }
                        const double weight = wy.weights[ny] * wz.weights[nz];
#pragma omp atomic update
                        fieldJ(cell_idx, j_idx, k_idx, X) +=
                            weight * charge * mpw * path / dt;
                    }
                }
                return std::make_pair(Bz * path, By * path);
            };

            double total_Bz = 0.0, total_By = 0.0;

            if (cell_diff == 0) {
                double path = end_x - start_x;
                auto [Bz, By] = process_cell(initial_cell, path);
                total_Bz = Bz;
                total_By = By;
            } else {
                const double boundary_x =
                    xCellSize * (initial_cell + (cell_diff > 0) - GHOST_CELLS);
                const double path1 = boundary_x - start_x;
                const double path2 = end_x - boundary_x;

                auto [Bz1, By1] = process_cell(initial_cell, path1);
                auto [Bz2, By2] = process_cell(final_cell, path2);

                total_Bz = Bz1 + Bz2;
                total_By = By1 + By2;
            }

            double3 v1, v2;
            v1.y() = - charge / mass * total_Bz;
            v1.z() = charge / mass * total_By;




            const double initial_cell_coord =
                (static_cast<int>((start_x / xCellSize))) * xCellSize;
            const double final_cell_coord =
                (static_cast<int>((end_x / xCellSize))) * xCellSize;

            if (final_cell == initial_cell) {
                double dv = charge / mass * (end_x - start_x);
                const int j_idx = wy.indices[0];
                const int k_idx = wz.indices[0];
                double Bz =
                    (fieldB(initial_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(initial_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By =
                    (fieldB(initial_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(initial_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                v2.y() = - dv * Bz;
                v2.z() = dv * By;
                // fieldJ(initial_cell, j_idx, k_idx, X) +=
                //     wy.weights[0] * wz.weights[0] * charge * mpw *
                //     particle.velocity.x();
                // fieldJ(initial_cell, j_idx + 1, k_idx, X) +=
                //     wy.weights[1] * wz.weights[0] * charge * mpw *
                //     particle.velocity.x();
                // fieldJ(initial_cell, j_idx, k_idx + 1, X) +=
                //     wy.weights[0] * wz.weights[1] * charge * mpw *
                //     particle.velocity.x();
                // fieldJ(initial_cell, j_idx + 1, k_idx + 1, X) +=
                //     wy.weights[1] * wz.weights[1] * charge * mpw *
                //     particle.velocity.x();
            }
            if (final_cell < initial_cell) {
                double dv = charge / mass;
                const int j_idx = wy.indices[0];
                const int k_idx = wz.indices[0];
                double Bz1 =
                    (fieldB(initial_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(initial_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By1 =
                    (fieldB(initial_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(initial_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                double Bz2 =
                    (fieldB(final_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(final_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By2 =
                    (fieldB(final_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(final_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                v2.y() = -(
                    dv * Bz1 * (initial_cell_coord - start_x) +
                    dv * Bz2 * (end_x - initial_cell_coord));

                v2.z() += dv * By1 * (initial_cell_coord - start_x) +
                                         dv * By2 * (end_x - initial_cell_coord);
                double dv2 = (end_x - initial_cell_coord) / dt;
                double dv1 = (initial_cell_coord - start_x) / dt;
                // fieldJ(initial_cell, j_idx, k_idx, X) +=
                //     wy.weights[0] * wz.weights[0] * charge * mpw * dv1;
                // fieldJ(initial_cell, j_idx + 1, k_idx, X) +=
                //     wy.weights[1] * wz.weights[0] * charge * mpw * dv1;
                // fieldJ(initial_cell, j_idx, k_idx + 1, X) +=
                //     wy.weights[0] * wz.weights[1] * charge * mpw * dv1;
                // fieldJ(initial_cell, j_idx + 1, k_idx + 1, X) +=
                //     wy.weights[1] * wz.weights[1] * charge * mpw * dv1;
                // fieldJ(final_cell, j_idx, k_idx, X) +=
                //     wy.weights[0] * wz.weights[0] * charge * mpw * dv2;
                // fieldJ(final_cell, j_idx + 1, k_idx, X) +=
                //     wy.weights[1] * wz.weights[0] * charge * mpw * dv2;
                // fieldJ(final_cell, j_idx, k_idx + 1, X) +=
                //     wy.weights[0] * wz.weights[1] * charge * mpw * dv2;
                // fieldJ(final_cell, j_idx + 1, k_idx + 1, X) +=
                //     wy.weights[1] * wz.weights[1] * charge * mpw * dv2;
            }
            if (final_cell > initial_cell) {
                double dv = charge / mass;
                const int j_idx = wy.indices[0];
                const int k_idx = wz.indices[0];
                double Bz1 =
                    (fieldB(initial_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(initial_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By1 =
                    (fieldB(initial_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(initial_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                double Bz2 =
                    (fieldB(final_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(final_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By2 =
                    (fieldB(final_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(final_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                v2.y() = -
                    (dv * Bz1 * (final_cell_coord - start_x) +
                    dv * Bz2 * (end_x - final_cell_coord));

                v2.z() =
                    dv * By1 * (final_cell_coord - start_x) +
                    dv * By2 * (end_x - final_cell_coord);
                // double dv2 = (end_x - final_cell_coord) / dt;
                // double dv1 = (final_cell_coord - start_x) / dt;
                // fieldJ(initial_cell, j_idx, k_idx, X) +=
                //     wy.weights[0] * wz.weights[0] * charge * mpw * dv1;
                // fieldJ(initial_cell, j_idx + 1, k_idx, X) +=
                //     wy.weights[1] * wz.weights[0] * charge * mpw * dv1;
                // fieldJ(initial_cell, j_idx, k_idx + 1, X) +=
                //     wy.weights[0] * wz.weights[1] * charge * mpw * dv1;
                // fieldJ(initial_cell, j_idx + 1, k_idx + 1, X) +=
                //     wy.weights[1] * wz.weights[1] * charge * mpw * dv1;
                // fieldJ(final_cell, j_idx, k_idx, X) +=
                //     wy.weights[0] * wz.weights[0] * charge * mpw * dv2;
                // fieldJ(final_cell, j_idx + 1, k_idx, X) +=
                //     wy.weights[1] * wz.weights[0] * charge * mpw * dv2;
                // fieldJ(final_cell, j_idx, k_idx + 1, X) +=
                //     wy.weights[0] * wz.weights[1] * charge * mpw * dv2;
                // fieldJ(final_cell, j_idx + 1, k_idx + 1, X) +=
                //     wy.weights[1] * wz.weights[1] * charge * mpw * dv2;
            }
            if(fabs(v1.y() - v2.y()) > 1e-16 || fabs(v1.z()- v2.z()) > 1e-16) {
                std::cout << start_x << " " << end_x << " " << initial_cell << " " << final_cell << " " << initial_cell_coord << " " << final_cell_coord << std::endl;
                std::cout << "v1.y() != v2.y() || v1.z() != v2.z() " <<  v1 << " " << v2 << std::endl;
            }
            particle.velocity.y() -= charge / mass * total_Bz;
            particle.velocity.z() += charge / mass * total_By;
            particle.coord.x() = end_x;
      //      std::cout << "x: " << particle.velocity << "\n";
        }
    }
}

void SimulationSymplectic::move_y(std::unique_ptr<ParticlesArray>& species,
                                  const Field3d& fieldB, Field3d& fieldJ,
                                  const double dt) {
    if (species->is_neutral())
        return;
    double xCellSize = domain.cell_size().x();
    double yCellSize = domain.cell_size().y();
    double zCellSize = domain.cell_size().z();
    double mass = species->mass();
    double charge = species->charge;
    double mpw = species->mpw();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            const double start_y = particle.coord.y();
            const double end_y = start_y + particle.velocity.y() * dt;
            const int initial_cell =
                static_cast<int>((start_y / yCellSize) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_y / yCellSize) + GHOST_CELLS);
            const int cell_diff = final_cell - initial_cell;

            const auto wx = compute_weights(
                (particle.coord.x() / xCellSize) + GHOST_CELLS, 0.0);
            const auto wz = compute_weights(
                (particle.coord.z() / zCellSize) + GHOST_CELLS, 0.0);

            const auto process_cell = [&](int cell_idx, double path) {
                const int i_idx = wx.indices[0];
                const int k_idx = wz.indices[0];
                double Bz =
                    (fieldB(i_idx, cell_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(i_idx, cell_idx, k_idx + 1, Z) * wz.weights[1]);
                double Bx =
                    (fieldB(i_idx, cell_idx, k_idx, X) * wx.weights[0] +
                     fieldB(i_idx + 1, cell_idx, k_idx, X) * wx.weights[1]);

#pragma omp simd collapse(2)
                for (int nx = 0; nx < 2; ++nx) {
                    for (int nz = 0; nz < 2; ++nz) {
                        const int i_idx = wx.indices[nx];
                        const int k_idx = wz.indices[nz];
                        // if (cell_idx < 0 || cell_idx >= 23 || i_idx < 0 ||
                        //     i_idx >= 23 || k_idx < 0 || k_idx >= 23) {
                        //     std::cout
                        //         << "Particle out of bounds: " << final_cell
                        //         << std::endl;
                        //     std::cout << particle << std::endl;
                        // }
                        const double weight = wx.weights[nx] * wz.weights[nz];
#pragma omp atomic update
                        fieldJ(i_idx, cell_idx, k_idx, Y) +=
                            weight * charge * mpw * path / dt;
                    }
                }
                return std::make_pair(Bz * path, Bx * path);
            };

            double total_Bz = 0.0, total_Bx = 0.0;

            if (cell_diff == 0) {
                double path = end_y - start_y;
                auto [Bz, Bx] = process_cell(initial_cell, path);
                total_Bz = Bz;
                total_Bx = Bx;
            } else {
                const double boundary_y =
                    yCellSize * (initial_cell + (cell_diff > 0) - GHOST_CELLS);
                const double path1 = boundary_y - start_y;
                const double path2 = end_y - boundary_y;

                auto [Bz1, Bx1] = process_cell(initial_cell, path1);
                auto [Bz2, Bx2] = process_cell(final_cell, path2);

                total_Bz = Bz1 + Bz2;
                total_Bx = Bx1 + Bx2;
            }

            particle.velocity.z() -= charge / mass * total_Bx;
            particle.velocity.x() += charge / mass * total_Bz;
            particle.coord.y() = end_y;
        //    std::cout << "y: " << particle.velocity << "\n";
        }
    }
}

void SimulationSymplectic::move_z(std::unique_ptr<ParticlesArray>& species,
                                  const Field3d& fieldB, Field3d& fieldJ,
                                  const double dt) {
    if (species->is_neutral())
        return;
    double xCellSize = domain.cell_size().x();
    double yCellSize = domain.cell_size().y();
    double zCellSize = domain.cell_size().z();
    double mass = species->mass();
    double charge = species->charge;
    double mpw = species->mpw();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            const double start_z = particle.coord.z();
            const double end_z = start_z + particle.velocity.z() * dt;
            const int initial_cell =
                static_cast<int>((start_z / zCellSize) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_z / zCellSize) + GHOST_CELLS);
            const int cell_diff = final_cell - initial_cell;

            const auto wx = compute_weights(
                (particle.coord.x() / xCellSize) + GHOST_CELLS, 0.0);
            const auto wy = compute_weights(
                (particle.coord.y() / yCellSize) + GHOST_CELLS, 0.0);

            const auto process_cell = [&](int cell_idx, double path) {
                const int i_idx = wx.indices[0];
                const int j_idx = wy.indices[0];
                double By = (fieldB(i_idx, j_idx, cell_idx, Y) * wy.weights[0] +
                             fieldB(i_idx, j_idx + 1, cell_idx, Y) *
                                 wy.weights[1]);
                double Bx =
                    (fieldB(i_idx, j_idx, cell_idx, X) * wx.weights[0] +
                     fieldB(i_idx + 1, j_idx, cell_idx, X) * wx.weights[1]);

#pragma omp simd collapse(2)
                for (int nx = 0; nx < 2; ++nx) {
                    for (int ny = 0; ny < 2; ++ny) {
                        const int i_idx = wx.indices[nx];
                        const int j_idx = wy.indices[ny];
                        // if (cell_idx < 0 || cell_idx >= 23 || j_idx < 0 ||
                        //     j_idx >= 23 || i_idx < 0 || i_idx >= 23) {
                        //     std::cout
                        //         << "Particle out of bounds: " << final_cell
                        //         << std::endl;
                        //     std::cout << particle << std::endl;
                        // }
                        const double weight = wx.weights[nx] * wy.weights[ny];
#pragma omp atomic update
                        fieldJ(i_idx, j_idx, cell_idx, Z) +=
                            weight * charge * mpw * path / dt;
                    }
                }
                return std::make_pair(By * path, Bx * path);
            };

            double total_By = 0.0, total_Bx = 0.0;

            if (cell_diff == 0) {
                double path = end_z - start_z;
                auto [By, Bx] = process_cell(initial_cell, path);
                total_By = By;
                total_Bx = Bx;
            } else {
                const double boundary_z =
                    zCellSize * (initial_cell + (cell_diff > 0) - GHOST_CELLS);
                const double path1 = boundary_z - start_z;
                const double path2 = end_z - boundary_z;

                auto [By1, Bx1] = process_cell(initial_cell, path1);
                auto [By2, Bx2] = process_cell(final_cell, path2);

                total_By = By1 + By2;
                total_Bx = Bx1 + Bx2;
            }

            particle.velocity.x() -= charge / mass * total_By;
            particle.velocity.y() += charge / mass * total_Bx;
            particle.coord.z() = end_z;
         //   std::cout << "z: " << particle.velocity << "\n";
        }
    }
}

void SimulationSymplectic::update_velocity(
    std::unique_ptr<ParticlesArray>& species, const double dt) {
    if (species->is_neutral())
        return;
    double xCellSize = domain.cell_size().x();
    double yCellSize = domain.cell_size().y();
    double zCellSize = domain.cell_size().z();
    double mass = species->mass();
    double charge = species->charge;
    double mpw = species->mpw();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            const auto wx = compute_weights(
                (particle.coord.x() / xCellSize) + GHOST_CELLS, 0.0);
            const auto wy = compute_weights(
                (particle.coord.y() / yCellSize) + GHOST_CELLS, 0.0);
            const auto wz = compute_weights(
                (particle.coord.z() / zCellSize) + GHOST_CELLS, 0.0);
            double Ex = 0.0, Ey = 0.0, Ez = 0.0;
#pragma omp simd collapse(2)
            for (int nx = 0; nx < 2; ++nx) {
                for (int ny = 0; ny < 2; ++ny) {
                    const int i_idx = wx.indices[nx];
                    const int j_idx = wy.indices[ny];
                    const int k_idx = wz.indices[0];
                    const double weight = wx.weights[nx] * wy.weights[ny];
                    Ez += fieldE(i_idx, j_idx, k_idx, Z) * weight;
                }
            }
#pragma omp simd collapse(2)
            for (int nx = 0; nx < 2; ++nx) {
                for (int nz = 0; nz < 2; ++nz) {
                    const int i_idx = wx.indices[nx];
                    const int k_idx = wz.indices[nz];
                    const int j_idx = wy.indices[0];
                    const double weight = wx.weights[nx] * wz.weights[nz];
                    Ey += fieldE(i_idx, j_idx, k_idx, Y) * weight;
                }
            }
#pragma omp simd collapse(2)
            for (int ny = 0; ny < 2; ++ny) {
                for (int nz = 0; nz < 2; ++nz) {
                    const int j_idx = wy.indices[ny];
                    const int k_idx = wz.indices[nz];
                    const int i_idx = wx.indices[0];
                    const double weight = wy.weights[ny] * wz.weights[nz];
                    Ex += fieldE(i_idx, j_idx, k_idx, X) * weight;
                }
            }
            particle.velocity.x() += Ex * charge * dt / mass;
            particle.velocity.y() += Ey * charge * dt / mass;
            particle.velocity.z() += Ez * charge * dt / mass;
        }
    }
}

void SimulationSymplectic::move_x_orig(std::unique_ptr<ParticlesArray>& species,
                                  const Field3d& fieldB, Field3d& fieldJ,
                                  const double dt) {
    if (species->is_neutral())
        return;
    double xCellSize = domain.cell_size().x();
    double yCellSize = domain.cell_size().y();
    double zCellSize = domain.cell_size().z();
    double mass = species->mass();
    double charge = species->charge;
    double mpw = species->mpw();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            const double start_x = particle.coord.x();
            const double end_x = start_x + particle.velocity.x() * dt;
            const int initial_cell =
                static_cast<int>((start_x / xCellSize) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_x / xCellSize) + GHOST_CELLS);
            const int initial_cell_coord =
                static_cast<int>((start_x / xCellSize)) * xCellSize;
            const int final_cell_coord =
                static_cast<int>((end_x / xCellSize)) * xCellSize;
            const int cell_diff = final_cell - initial_cell;
            
            const auto wy = compute_weights(
                (particle.coord.y() / yCellSize) + GHOST_CELLS, 0.0);
            const auto wz = compute_weights(
                (particle.coord.z() / zCellSize) + GHOST_CELLS, 0.0);

            if(final_cell == initial_cell) {
                double dv =  charge / mass * (end_x - start_x);
                const int j_idx = wy.indices[0];
                const int k_idx = wz.indices[0];
                double Bz =
                    (fieldB(initial_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(initial_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By =
                    (fieldB(initial_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(initial_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                particle.velocity.y() -= dv * Bz;
                particle.velocity.z() += dv * By;
                fieldJ(initial_cell, j_idx, k_idx, X) +=
                    wy.weights[0] * wz.weights[0] * charge * mpw *
                    particle.velocity.x();
                fieldJ(initial_cell, j_idx+1, k_idx, X) +=
                    wy.weights[1] * wz.weights[0] * charge * mpw *
                    particle.velocity.x();
                fieldJ(initial_cell, j_idx, k_idx+1, X) +=
                    wy.weights[0] * wz.weights[1] * charge * mpw *
                    particle.velocity.x();
                fieldJ(initial_cell, j_idx+1, k_idx+1, X) +=
                    wy.weights[1] * wz.weights[1] * charge * mpw *
                    particle.velocity.x();
            }
            if(final_cell < initial_cell) {
                double dv = charge / mass;
                const int j_idx = wy.indices[0];
                const int k_idx = wz.indices[0];
                double Bz1 =
                    (fieldB(initial_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(initial_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By1 =
                    (fieldB(initial_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(initial_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                double Bz2 =
                    (fieldB(final_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(final_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By2 =
                    (fieldB(final_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(final_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                particle.velocity.y() -=
                    dv * Bz1 * (initial_cell_coord - start_x) +
                    dv * Bz2 * (end_x - initial_cell_coord);

                particle.velocity.z() += dv * By1 * (initial_cell_coord - start_x) +
                                         dv * By2 * (end_x - initial_cell_coord);
                double dv2 = (end_x - initial_cell_coord) / dt;
                double dv1 = (initial_cell_coord - start_x) / dt;
                fieldJ(initial_cell, j_idx, k_idx, X) +=
                    wy.weights[0] * wz.weights[0] * charge * mpw * dv1;
                fieldJ(initial_cell, j_idx + 1, k_idx, X) +=
                    wy.weights[1] * wz.weights[0] * charge * mpw * dv1;
                fieldJ(initial_cell, j_idx, k_idx + 1, X) +=
                    wy.weights[0] * wz.weights[1] * charge * mpw * dv1;
                fieldJ(initial_cell, j_idx + 1, k_idx + 1, X) +=
                    wy.weights[1] * wz.weights[1] * charge * mpw * dv1;
                fieldJ(final_cell, j_idx, k_idx, X) +=
                    wy.weights[0] * wz.weights[0] * charge * mpw * dv2;
                fieldJ(final_cell, j_idx + 1, k_idx, X) +=
                    wy.weights[1] * wz.weights[0] * charge * mpw * dv2;
                fieldJ(final_cell, j_idx, k_idx + 1, X) +=
                    wy.weights[0] * wz.weights[1] * charge * mpw * dv2;
                fieldJ(final_cell, j_idx + 1, k_idx + 1, X) +=
                    wy.weights[1] * wz.weights[1] * charge * mpw * dv2;
            }
            if (final_cell > initial_cell) {
                double dv = charge / mass;
                const int j_idx = wy.indices[0];
                const int k_idx = wz.indices[0];
                double Bz1 =
                    (fieldB(initial_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(initial_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By1 =
                    (fieldB(initial_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(initial_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                double Bz2 =
                    (fieldB(final_cell, j_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(final_cell, j_idx, k_idx + 1, Z) * wz.weights[1]);
                double By2 =
                    (fieldB(final_cell, j_idx, k_idx, Y) * wy.weights[0] +
                     fieldB(final_cell, j_idx + 1, k_idx, Y) * wy.weights[1]);

                particle.velocity.y() -=
                    dv * Bz1 * (final_cell_coord - start_x) +
                    dv * Bz2 * (end_x - final_cell_coord);

                particle.velocity.z() +=
                    dv * By1 * (final_cell_coord - start_x) +
                    dv * By2 * (end_x - final_cell_coord);
                double dv2 = (end_x - final_cell_coord) / dt;
                double dv1 = (final_cell_coord - start_x) / dt;
                fieldJ(initial_cell, j_idx, k_idx, X) +=
                    wy.weights[0] * wz.weights[0] * charge * mpw * dv1;
                fieldJ(initial_cell, j_idx + 1, k_idx, X) +=
                    wy.weights[1] * wz.weights[0] * charge * mpw * dv1;
                fieldJ(initial_cell, j_idx, k_idx + 1, X) +=
                    wy.weights[0] * wz.weights[1] * charge * mpw * dv1;
                fieldJ(initial_cell, j_idx + 1, k_idx + 1, X) +=
                    wy.weights[1] * wz.weights[1] * charge * mpw * dv1;
                fieldJ(final_cell, j_idx, k_idx, X) +=
                    wy.weights[0] * wz.weights[0] * charge * mpw * dv2;
                fieldJ(final_cell, j_idx + 1, k_idx, X) +=
                    wy.weights[1] * wz.weights[0] * charge * mpw * dv2;
                fieldJ(final_cell, j_idx, k_idx + 1, X) +=
                    wy.weights[0] * wz.weights[1] * charge * mpw * dv2;
                fieldJ(final_cell, j_idx + 1, k_idx + 1, X) +=
                    wy.weights[1] * wz.weights[1] * charge * mpw * dv2;
            }

            particle.coord.x() = end_x;
        }
    }
}

void SimulationSymplectic::move_y_orig(std::unique_ptr<ParticlesArray>& species,
                                       const Field3d& fieldB, Field3d& fieldJ,
                                       const double dt) {
    if (species->is_neutral())
        return;
    double xCellSize = domain.cell_size().x();
    double yCellSize = domain.cell_size().y();
    double zCellSize = domain.cell_size().z();
    double mass = species->mass();
    double charge = species->charge;
    double mpw = species->mpw();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            const double start_y = particle.coord.y();
            const double end_y = start_y + particle.velocity.y() * dt;
            const int initial_cell =
                static_cast<int>((start_y / yCellSize) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_y / yCellSize) + GHOST_CELLS);
            const int cell_diff = final_cell - initial_cell;

            const auto wx = compute_weights(
                (particle.coord.x() / xCellSize) + GHOST_CELLS, 0.0);
            const auto wz = compute_weights(
                (particle.coord.z() / zCellSize) + GHOST_CELLS, 0.0);

            const auto process_cell = [&](int cell_idx, double path) {
                const int i_idx = wx.indices[0];
                const int k_idx = wz.indices[0];
                double Bz =
                    (fieldB(i_idx, cell_idx, k_idx, Z) * wz.weights[0] +
                     fieldB(i_idx, cell_idx, k_idx + 1, Z) * wz.weights[1]);
                double Bx =
                    (fieldB(i_idx, cell_idx, k_idx, X) * wx.weights[0] +
                     fieldB(i_idx + 1, cell_idx, k_idx, X) * wx.weights[1]);

#pragma omp simd collapse(2)
                for (int nx = 0; nx < 2; ++nx) {
                    for (int nz = 0; nz < 2; ++nz) {
                        const int i_idx = wx.indices[nx];
                        const int k_idx = wz.indices[nz];
                        if (cell_idx < 0 || cell_idx >= 23 || i_idx < 0 ||
                            i_idx >= 23 || k_idx < 0 || k_idx >= 23) {
                            std::cout
                                << "Particle out of bounds: " << final_cell
                                << std::endl;
                            std::cout << particle << std::endl;
                        }
                        const double weight = wx.weights[nx] * wz.weights[nz];
#pragma omp atomic update
                        fieldJ(i_idx, cell_idx, k_idx, Y) +=
                            weight * charge * mpw * path / dt;
                    }
                }
                return std::make_pair(Bz * path, Bx * path);
            };

            double total_Bz = 0.0, total_Bx = 0.0;

            if (cell_diff == 0) {
                double path = end_y - start_y;
                auto [Bz, Bx] = process_cell(initial_cell, path);
                total_Bz = Bz;
                total_Bx = Bx;
            } else {
                const double boundary_y =
                    yCellSize * (initial_cell + (cell_diff > 0) - GHOST_CELLS);
                const double path1 = boundary_y - start_y;
                const double path2 = end_y - boundary_y;

                auto [Bz1, Bx1] = process_cell(initial_cell, path1);
                auto [Bz2, Bx2] = process_cell(final_cell, path2);

                total_Bz = Bz1 + Bz2;
                total_Bx = Bx1 + Bx2;
            }

            particle.velocity.z() -= charge / mass * total_Bx;
            particle.velocity.x() += charge / mass * total_Bz;
            particle.coord.y() = end_y;
        }
    }
}

void SimulationSymplectic::move_z_orig(std::unique_ptr<ParticlesArray>& species,
                                       const Field3d& fieldB, Field3d& fieldJ,
                                       const double dt) {
    if (species->is_neutral())
        return;
    double xCellSize = domain.cell_size().x();
    double yCellSize = domain.cell_size().y();
    double zCellSize = domain.cell_size().z();
    double mass = species->mass();
    double charge = species->charge;
    double mpw = species->mpw();

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < species->size(); ++pk) {
        for (auto& particle : species->particlesData(pk)) {
            const double start_z = particle.coord.z();
            const double end_z = start_z + particle.velocity.z() * dt;
            const int initial_cell =
                static_cast<int>((start_z / zCellSize) + GHOST_CELLS);
            const int final_cell =
                static_cast<int>((end_z / zCellSize) + GHOST_CELLS);
            const int cell_diff = final_cell - initial_cell;

            const auto wx = compute_weights(
                (particle.coord.x() / xCellSize) + GHOST_CELLS, 0.0);
            const auto wy = compute_weights(
                (particle.coord.y() / yCellSize) + GHOST_CELLS, 0.0);

            const auto process_cell = [&](int cell_idx, double path) {
                const int i_idx = wx.indices[0];
                const int j_idx = wy.indices[0];
                double By =
                    (fieldB(i_idx, j_idx, cell_idx, Y) * wy.weights[0] +
                     fieldB(i_idx, j_idx + 1, cell_idx, Y) * wy.weights[1]);
                double Bx =
                    (fieldB(i_idx, j_idx, cell_idx, X) * wx.weights[0] +
                     fieldB(i_idx + 1, j_idx, cell_idx, X) * wx.weights[1]);

#pragma omp simd collapse(2)
                for (int nx = 0; nx < 2; ++nx) {
                    for (int ny = 0; ny < 2; ++ny) {
                        const int i_idx = wx.indices[nx];
                        const int j_idx = wy.indices[ny];
                        if (cell_idx < 0 || cell_idx >= 23 || j_idx < 0 ||
                            j_idx >= 23 || i_idx < 0 || i_idx >= 23) {
                            std::cout
                                << "Particle out of bounds: " << final_cell
                                << std::endl;
                            std::cout << particle << std::endl;
                        }
                        const double weight = wx.weights[nx] * wy.weights[ny];
#pragma omp atomic update
                        fieldJ(i_idx, j_idx, cell_idx, Z) +=
                            weight * charge * mpw * path / dt;
                    }
                }
                return std::make_pair(By * path, Bx * path);
            };

            double total_By = 0.0, total_Bx = 0.0;

            if (cell_diff == 0) {
                double path = end_z - start_z;
                auto [By, Bx] = process_cell(initial_cell, path);
                total_By = By;
                total_Bx = Bx;
            } else {
                const double boundary_z =
                    zCellSize * (initial_cell + (cell_diff > 0) - GHOST_CELLS);
                const double path1 = boundary_z - start_z;
                const double path2 = end_z - boundary_z;

                auto [By1, Bx1] = process_cell(initial_cell, path1);
                auto [By2, Bx2] = process_cell(final_cell, path2);

                total_By = By1 + By2;
                total_Bx = Bx1 + Bx2;
            }

            particle.velocity.x() -= charge / mass * total_By;
            particle.velocity.y() += charge / mass * total_Bx;
            particle.coord.z() = end_z;
        }
    }
}
