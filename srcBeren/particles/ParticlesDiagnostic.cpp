#include <functional>
#include <iterator>

#include "ParticlesArray.h"
#include "collision.h"
#include "containers.h"
#include "coupling/shape.h"
#include "grid/domain.h"
#include "sgs.h"
#include "timer.h"

std::ostream& operator<<(std::ostream& out, const ParticleSimple& particle) {
    out << particle.coord << " " << particle.velocity;
    return out;
}
std::ostream& operator<<(std::ostream& out, const ParticleMass& particle) {
    out << particle.coord << " " << particle.velocity << " " << particle.mass;
    return out;
}

double PulseFromKev(double kev, double mass) {
    double gama = kev / SGS::MC2 + mass;
    return sqrt((gama * gama) - mass * mass);
}

ParticleMoments ParticlesArray::diagnostic_moments() const {
    RECORD_TIMER;
    const Vector3I cells = particlesData.size();
    Array2D<ParticleMoments> rowMoments(cells.x(), cells.y());
    const double factor = 0.5 * mass_ * mpw_;
#pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < cells.x(); ++i) {
        for (int j = 0; j < cells.y(); ++j) {
            double ex = 0, ey = 0, ez = 0;
            int64_t count = 0;
            for (int k = 0; k < cells.z(); ++k) {
                count += std::ssize(particlesData(i, j, k));
                for (const Particle& particle : particlesData(i, j, k)) {
                    const Vector3R& v = particle.velocity;
                    ex += factor * v.x() * v.x();
                    ey += factor * v.y() * v.y();
                    ez += factor * v.z() * v.z();
                }
            }
            rowMoments(i, j) = {Vector3R(ex, ey, ez), count};
        }
    }
    // Sum rows in mesh order: reproducible for any thread count, with O(Nx*Ny) storage.
    double ex = 0, ey = 0, ez = 0;
    int64_t count = 0;
    for (int i = 0; i < cells.x(); ++i) {
        for (int j = 0; j < cells.y(); ++j) {
            const ParticleMoments& m = rowMoments(i, j);
            ex += m.kinetic.x();
            ey += m.kinetic.y();
            ez += m.kinetic.z();
            count += m.count;
        }
    }
    return {Vector3R(ex, ey, ez), count};
}

double ParticlesArray::get_kinetic_energy() const {
    RECORD_TIMER;
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            energy += get_energy_particle(particle.velocity, mass_, mpw_);
        }
    }

    return energy;
}

double ParticlesArray::get_kinetic_energy(int dim) const {
    RECORD_TIMER;
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            Vector3R velocity(0, 0, 0);
            velocity[dim] = particle.velocity[dim];
            energy += get_energy_particle(velocity, mass_, mpw_);
        }
    }

    return energy;
}

Vector3R ParticlesArray::get_kinetic_energy_component() const {
    RECORD_TIMER;
    double enx = 0;
    double eny = 0;
    double enz = 0;
#pragma omp parallel for reduction(+ : enx) reduction(+ : eny) reduction(+ : enz)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            auto velocity = particle.velocity;
            enx += get_energy_particle(velocity.x(), mass_, mpw_);
            eny += get_energy_particle(velocity.y(), mass_, mpw_);
            enz += get_energy_particle(velocity.z(), mass_, mpw_);
        }
    }

    return Vector3R(enx, eny, enz);
}
double ParticlesArray::get_kinetic_energy(int dim1, int dim2) const {
    RECORD_TIMER;
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            Vector3R velocity(0, 0, 0);
            velocity[dim1] = particle.velocity[dim1];
            velocity[dim2] = particle.velocity[dim2];
            energy += get_energy_particle(velocity, mass_, mpw_);
        }
    }

    return energy;
}

EnergySpectrum ParticlesArray::calculate_energy_spectrum() const {
    RECORD_TIMER;
    const int num_bins = 1000;
    // Находим минимальную и максимальную энергии
    double min_energy = std::numeric_limits<double>::max();
    double max_energy = std::numeric_limits<double>::lowest();

#pragma omp parallel for reduction(min : min_energy) reduction(max : max_energy)
    for (int k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            double e = get_energy_particle(particle.velocity, mass_, mpw_);
            if (e < min_energy)
                min_energy = e;
            if (e > max_energy)
                max_energy = e;
        }
    }

    // Обработка случая с одинаковыми энергиями
    if (min_energy == max_energy) {
        min_energy -= 0.1;
        max_energy += 0.1;
    }
    double bin_width = (max_energy - min_energy) / num_bins;

    // Счётчики частиц в бинах
    std::vector<int> counts(num_bins, 0);

// Заполнение гистограммы
#pragma omp parallel
    {
        std::vector<int> local_counts(num_bins, 0);

#pragma omp for
        for (int k = 0; k < size(); ++k) {
            for (const auto& particle : particlesData(k)) {
                double e = get_energy_particle(particle.velocity, mass_, mpw_);
                int index = static_cast<int>((e - min_energy) / bin_width);
                index = std::clamp(index, 0, num_bins - 1);
                local_counts[index]++;
            }
        }

#pragma omp critical
        {
            for (int i = 0; i < num_bins; ++i) {
                counts[i] += local_counts[i];
            }
        }
    }
    return EnergySpectrum(counts, min_energy, max_energy);
}
