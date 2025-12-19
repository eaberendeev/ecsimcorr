#include "interpolation.h"
#include "simulation_ecsim_corr.h"

void SimulationEcsimCorr::correctv(ParticlesArray& sort, const double dt) {
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
            const Vector3R end = particle.coord;
            const auto initVelocity = particle.initVelocity;
            const auto velocity = particle.velocity;
            const Vector3R coord = end - 0.5 * dt * velocity;
            const auto norm_coord = sort.to_cell_coordinates(coord);
            const Vector3R Ep = interpolateE(fieldEp, norm_coord, SHAPE);

            const Vector3R v12 = 0.5 * (velocity + initVelocity);

            jp_cell_loc += mpw * charge * v12.dot(Ep);
        }
        jp_cell += jp_cell_loc;
    }

    const double energyJeEn = calc_JE(fieldEn, currentOnGrid, bounds);
    const double energyJeE = calc_JE(fieldE, currentOnGrid, bounds);

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

    std::cout << "lambda " << lambda << " " << lambda * lambda << "\n";
}
