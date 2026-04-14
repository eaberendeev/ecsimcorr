#include "interpolation.h"
#include "simulation_ecsim_corr.h"

void SimulationEcsimCorr::correctv(ParticlesArray& sort, const double dt) {
    if (sort.is_neutral())
        return;

    const Field3d fieldEp_full = fieldEp + fieldE_external;
    const Field3d fieldEp_corr_full = 0.5*(fieldE + fieldEn) + fieldE_external;

    const double charge = sort.charge;
    const double mpw = sort.mpw();
    const auto& currentOnGrid = sort.currentOnGrid;
    const auto& domain = sort.get_domain();

    double jp_cell = 0;
#pragma omp parallel for schedule(guided) reduction(+ : jp_cell)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        const auto& particles = sort.particlesData(pk);
        double jp_cell_loc = 0;
        for (const auto& particle : particles) {
            const Vector3R end = particle.coord;
            const auto initVelocity = particle.initVelocity;
            const auto velocity = particle.velocity;
            const Vector3R coord = end - 0.5 * dt * velocity;
            const auto norm_coord = domain.to_cell_coordinates(coord);
            const Vector3R Ep = interpolateE(fieldEp_full, norm_coord, SHAPE);

            const Vector3R v12 = 0.5 * (velocity + initVelocity);

            jp_cell_loc += mpw * charge * v12.dot(Ep);
        }
        jp_cell += jp_cell_loc;
    }

    const double energyJe_corr = calc_JE(fieldEp_corr_full, currentOnGrid);

    // change to
    // energy += get_energy_particle(particle.velocity,
    // mass_, mpw_);

    const double energyK = sort.get_kinetic_energy();
    const double lambda = sqrt(1 + dt * (energyJe_corr - jp_cell) / energyK);

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        for (auto& particle : sort.particlesData(pk)) {
            particle.velocity = lambda * particle.velocity;
        }
    }

    std::cout << "lambda " << lambda << " " << lambda * lambda << "\n";
}

// using ShapeFunction = double (*)(const double&);

// template void move_and_calc_current_impl<Shape, 2>(ParticlesArray& particles,
//     const double dt, Field3d& fieldJ);
// template void move_and_calc_current_impl<Shape2, 2>(ParticlesArray& particles,
//                                                     const double dt,
//                                                     Field3d& fieldJ);

// void move_and_calc_current(ParticlesArray& particles,
//                                            const double dt, Field3d& fieldJ,
//                                            ShapeType type) {
//     if (is_neutral())
//         return;

//     switch (type) {
//         case ShapeType::NGP:
//             std::cout << "Move and calc current for NGP is not supported\n"
//                       << std::endl;
//             exit(-1);
//         case ShapeType::Linear:
//             move_and_calc_current_impl<Shape, 2>(particles, dt, fieldJ);
//             break;
//         case ShapeType::Quadratic:
//             move_and_calc_current_impl<Shape2, 2>(particles, dt, fieldJ);
//             break;
//     }
// }

// template <ParticlesArray::ShapeFunction ShapeFn, int ShapeSize>
// void move_and_calc_current_impl(ParticlesArray& particles, Domain& domain,
//                                                 const double dt,
//                                                 Field3d& fieldJ) {
//     constexpr auto SMAX = 2 * ShapeSize;

//     const double qx = charge * domain_.cell_size().x() / (6 * dt) * mpw_;
//     const double qy = charge * domain_.cell_size().y() / (6 * dt) * mpw_;
//     const double qz = charge * domain_.cell_size().z() / (6 * dt) * mpw_;

// // TODO: change base_ to cell index from ParticlesData
// #pragma omp parallel for schedule(dynamic, 64)
//     for (auto pk = 0; pk < size(); ++pk) {
//         auto& particles = particlesData(pk);
//         if (particles.empty()) {
//             continue;
//         }

//         ParticleShape<ShapeFn, SMAX> start_shape;
//         ParticleShape<ShapeFn, SMAX> end_shape;
//         CurrentBuffer<SMAX> curBuf, cellBuf;
//         cellBuf.zero();
//         curBuf.zero();
//         start_shape.fill_zero();

//         for (auto& particle : particles) {
//             Vector3R start = particle.coord;
//             start_shape.fill_from_normalized(domain_.to_cell_coordinates(start),
//                                              GHOST_CELLS);
//             particle.move(dt);

//             Vector3R end = particle.coord;
//             end_shape.fill_from_normalized(domain_.to_cell_coordinates(end),
//                                            start_shape.base_, GHOST_CELLS);
//             decompose_esirkepov_current(start_shape, end_shape, qx, qy, qz,
//                                         curBuf);

//             cellBuf += curBuf;
//         }
//         auto [start_x, start_y, start_z] = start_shape.start_.split();
//         flush_current_buffer(fieldJ, cellBuf, start_x, start_y, start_z);
//     }
// }
