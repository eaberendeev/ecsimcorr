#include "simulation_ecsim_corr.h"

void SimulationEcsimCorr::correctv(ParticlesArray& sort, const Field3d& Jfull,
                                   const double dt) {
    if (sort.is_neutral())
        return;

    const double charge = sort.charge;
    const double mpw = sort.mpw();
    const auto& currentOnGrid = sort.currentOnGrid;

    double jp_cell = 0;
#pragma omp parallel for schedule(guided) reduction(+ : jp_cell)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        const auto& particles = sort.particlesData(pk);
        double jp_cell_loc = 0;
        for (const auto& particle : particles) {
            const double3 end = particle.coord;
            const auto initVelocity = particle.initVelocity;
            const auto velocity = particle.velocity;
            const double3 coord = end - 0.5 * dt * velocity;
            const auto norm_coord = sort.normalize_coord(coord);
            const double3 Ep = interpolateE(fieldEp, norm_coord, SHAPE);
            const double3 E = interpolateE(fieldE, norm_coord, SHAPE);

            const double3 v12 = 0.5 * (velocity + initVelocity);

            jp_cell_loc += 0.5 * mpw * charge * dot(v12, (Ep + E));
        }
        jp_cell += jp_cell_loc;
    }

    const double energyJeEn = calc_JE(fieldEn, currentOnGrid, bounds);
    const double energyJeE = calc_JE(fieldE, currentOnGrid, bounds);
    const double energyJpEp = calc_JE(fieldEp, Jfull, bounds);
    const double energyJpE = calc_JE(fieldE, Jfull, bounds);
    // change to
    // energy += get_energy_particle(particle.velocity,
    // _mass, _mpw);

    const double energyK = sort.get_kinetic_energy();
    const double lambda =
        sqrt(1 + dt * (0.5 * (energyJeEn + energyJeE) - jp_cell) / energyK);

#pragma omp parallel for schedule(dynamic, 64)
    for (auto pk = 0; pk < sort.size(); ++pk) {
        for (auto& particle : sort.particlesData(pk)) {
            particle.velocity = lambda * particle.velocity;
        }
    }

    const double energyK2 = sort.get_kinetic_energy();
    std::cout << "lambda " << lambda << " " << lambda * lambda << " "
              << energyK2 - energyK << " "
              << 0.5 * dt * (energyJeEn + energyJeE - energyJpEp - energyJpE)
              << "\n";
}
