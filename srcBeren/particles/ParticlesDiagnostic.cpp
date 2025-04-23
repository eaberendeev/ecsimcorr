#include <functional>

#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "service.h"
#include "sgs.h"

std::ostream& operator<<(std::ostream& out,const double3& val){
	  out <<  val.x() << " " << val.y() << " " << val.z();
	  return out;
} 

std::ostream& operator<<(std::ostream& out, const ParticleSimple& particle){
	out << particle.coord << " " << particle.velocity;
	return out;
} 
std::ostream& operator<<(std::ostream& out, const ParticleMass& particle){
	out << particle.coord << " " << particle.velocity << " " << particle.mass;
	return out;
} 



double PulseFromKev(double kev, double mass){
  double gama = kev / SGS::MC2 + mass;
  return sqrt((gama*gama)- mass);
}

// Template specializations need to be explicitly instantiated in the cpp file
template void ParticlesArray::density_on_grid_update_impl<Shape, 2>();
template void ParticlesArray::density_on_grid_update_impl<Shape2, 2>();
// Implementation of the template function
template <ParticlesArray::ShapeFunction ShapeFn, int ShapeSize>
void ParticlesArray::density_on_grid_update_impl() {
    constexpr auto SMAX = 2 * ShapeSize;
    densityOnGrid.setZero();
    // std::cout << "density_on_grid_update_impl " << SMAX << std::endl;

#pragma omp parallel for schedule(dynamic, 64)
    for (auto j = 0; j < size(); ++j) {
        alignas(64) double sx[SMAX];
        alignas(64) double sy[SMAX];
        alignas(64) double sz[SMAX];

        for (const auto& particle : particlesData(j)) {
            // Vectorizable coordinate calculations
            const double ix = particle.coord.x() / xCellSize;
            const double iy = particle.coord.y() / yCellSize;
            const double iz = particle.coord.z() / zCellSize;

            const int xk = int(ix);
            const int yk = int(iy);
            const int zk = int(iz);

// Vectorizable shape calculations
#pragma omp simd
            for (int n = 0; n < SMAX; ++n) {
                sx[n] = ShapeFn(-ix + double(xk - GHOST_CELLS + n));
                sy[n] = ShapeFn(-iy + double(yk - GHOST_CELLS + n));
                sz[n] = ShapeFn(-iz + double(zk - GHOST_CELLS + n));
            }

            const double weight = is_neutral() ? _mpw : _mpw * charge;

// Density accumulation with loop unrolling
// #pragma unroll
            for (int n = 0; n < SMAX; ++n) {
                const int indx = xk + n;
                const double sxw = sx[n];

                for (int m = 0; m < SMAX; ++m) {
                    const int indy = yk + m;
                    const double sxyw = sxw * sy[m];

                    for (int k = 0; k < SMAX; ++k) {
                        const int indz = zk + k;
#pragma omp atomic update
                        densityOnGrid(indx, indy, indz, 0) +=
                            weight * sxyw * sz[k];
                    }
                }
            }
        }
    }
}

void ParticlesArray::density_on_grid_update_impl_ngp() {
    densityOnGrid.setZero();
#pragma omp parallel for schedule(dynamic, 64)
    for (auto j = 0; j < size(); ++j) {

        for (const auto& particle : particlesData(j)) {
            // Vectorizable coordinate calculations
            const double x = particle.coord.x() / xCellSize;
            const double y = particle.coord.y() / yCellSize;
            const double z = particle.coord.z() / zCellSize;

            const int indx = ngp(x);
            const int indy = ngp(y);
            const int indz = ngp(z);

            const double weight = is_neutral() ? _mpw : _mpw * charge;

#pragma omp atomic update
            densityOnGrid(indx, indy, indz, 0) += weight;
        }
    }
}


// Public interface that selects appropriate implementation
void ParticlesArray::density_on_grid_update(
    ShapeType type) {
    switch (type) {
        case ShapeType::NGP:
            density_on_grid_update_impl_ngp();
            break;
        case ShapeType::Linear:
            density_on_grid_update_impl<Shape, 2>();
            break;
        case ShapeType::Quadratic:
            density_on_grid_update_impl<Shape2, 2>();
            break;
    }
}

double ParticlesArray::get_kinetic_energy() const{
    double energy = 0;
#pragma omp parallel for reduction(+:energy)
    for(auto k = 0; k < size(); ++k){
        for(const auto& particle : particlesData(k)){
            energy += get_energy_particle(particle.velocity, _mass, _mpw);
        }
    }

    return energy;
}

double ParticlesArray::get_init_kinetic_energy() const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            energy += get_energy_particle(particle.initVelocity, _mass, _mpw);
        }
    }

    return energy;
}

double ParticlesArray::get_kinetic_energy(int dim) const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            double3 velocity(0,0,0);
            velocity(dim) = particle.velocity(dim);
            energy +=
                get_energy_particle(velocity, _mass, _mpw);
        }
    }

    return energy;
}

double3 ParticlesArray::get_kinetic_energy_component() const {
    double enx = 0;
    double eny = 0;
    double enz = 0;
#pragma omp parallel for reduction(+ : enx) reduction(+ : eny) reduction(+ : enz)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            auto velocity = particle.velocity;
            enx += get_energy_particle(velocity.x(), _mass, _mpw);
            eny += get_energy_particle(velocity.y(), _mass, _mpw);
            enz += get_energy_particle(velocity.z(), _mass, _mpw);
        }
    }

    return double3(enx, eny, enz);
}
double ParticlesArray::get_kinetic_energy(int dim1, int dim2) const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            double3 velocity(0, 0, 0);
            velocity(dim1) = particle.velocity(dim1);
            velocity(dim2) = particle.velocity(dim2);
            energy += get_energy_particle(velocity, _mass, _mpw);
        }
    }

    return energy;
}

EnergySpectrum ParticlesArray::calculate_energy_spectrum() const {
    const int num_bins = 1000;
    // Находим минимальную и максимальную энергии
    double min_energy = std::numeric_limits<double>::max();
    double max_energy = std::numeric_limits<double>::lowest();

#pragma omp parallel for reduction(min : min_energy) reduction(max : max_energy)
    for (int k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            double e = get_energy_particle(particle.velocity, _mass, _mpw);
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
                double e = get_energy_particle(particle.velocity, _mass, _mpw);
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
