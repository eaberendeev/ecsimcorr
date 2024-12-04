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
    densityOnGrid.set_zero();
    // std::cout << "density_on_grid_update_impl " << SMAX << std::endl;

#pragma omp parallel for
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

            const double weight = _mpw * charge;

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
    densityOnGrid.set_zero();
#pragma omp parallel for
    for (auto j = 0; j < size(); ++j) {

        for (const auto& particle : particlesData(j)) {
            // Vectorizable coordinate calculations
            const double x = particle.coord.x() / xCellSize;
            const double y = particle.coord.y() / yCellSize;
            const double z = particle.coord.z() / zCellSize;

            const int indx = ngp(x);
            const int indy = ngp(y);
            const int indz = ngp(z);

            const double weight = _mpw * charge;

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

void ParticlesArray::get_P() {
    Pxx.set_zero();
    Pyy.set_zero();
    Pzz.set_zero();

    constexpr auto SMAX = SHAPE_SIZE;

#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        double wx[2];
        double wy[2];
        double wz[2];
        for (auto& particle : particlesData(pk)) {
            const auto coord = particle.coord;
            const auto velocity = particle.velocity;

            auto x = coord.x() / xCellSize + GHOST_CELLS;
            auto y = coord.y() / yCellSize + GHOST_CELLS;
            auto z = coord.z() / zCellSize + GHOST_CELLS;

            const auto intx = int(x);
            const auto inty = int(y);
            const auto intz = int(z);

            wx[1] = (x - intx);
            wx[0] = 1 - wx[1];
            wy[1] = (y - inty);
            wy[0] = 1 - wy[1];
            wz[1] = (z - intz);
            wz[0] = 1 - wz[1];

            double vx = velocity.x();
            double vy = velocity.y();
            double vz = velocity.z();
            double pxx = _mass * vx * vx * _mpw;
            double pyy = _mass * vy * vy * _mpw;
            double pzz = _mass * vz * vz * _mpw;

            for (int nx = 0; nx < SMAX; ++nx) {
                const int i = intx + nx;
                for (int ny = 0; ny < SMAX; ++ny) {
                    const int j = inty + ny;
                    for (int nz = 0; nz < SMAX; ++nz) {
                        const int k = intz + nz;
                        const auto sx = wx[nx] * wy[ny] * wz[nz];
#pragma omp atomic update
                        Pxx(i, j, k, 0) += sx * pxx;
#pragma omp atomic update
                        Pyy(i, j, k, 0) += sx * pyy;
#pragma omp atomic update
                        Pzz(i, j, k, 0) += sx * pzz;
                    }
                }
            }
        }
    }
}

void ParticlesArray::get_Pr() {
    Pxx.set_zero();
    Pyy.set_zero();
    Pzz.set_zero();

    constexpr auto SMAX = SHAPE_SIZE;

#pragma omp parallel for
    for (auto pk = 0; pk < size(); ++pk) {
        double wx[2];
        double wy[2];
        double wz[2];
        for (auto& particle : particlesData(pk)) {
            const auto coord = particle.coord;
            const auto velocity = particle.velocity;

            auto x = coord.x() / xCellSize + GHOST_CELLS;
            auto y = coord.y() / yCellSize + GHOST_CELLS;
            auto z = coord.z() / zCellSize + GHOST_CELLS;

            const auto intx = int(x);
            const auto inty = int(y);
            const auto intz = int(z);

            wx[1] = (x - intx);
            wx[0] = 1 - wx[1];
            wy[1] = (y - inty);
            wy[0] = 1 - wy[1];
            wz[1] = (z - intz);
            wz[0] = 1 - wz[1];

            double x0 =
                0.5 * xCellSize * xCellCount;   // to do: domain_get_center
            double R = sqrt((coord.x() - x0) * (coord.x() - x0) +
                            (coord.y() - x0) * (coord.y() - x0));

            double vr = ((coord.x() - x0) / R) * velocity.x() +
                        ((coord.y() - x0) / R) * velocity.y();
            double vp = -((coord.y() - x0) / R) * velocity.x() +
                        ((coord.x() - x0) / R) * velocity.y();
            double vz = velocity.z();
            double prr = _mass * vr * vr * _mpw;
            double ppp = _mass * vp * vp * _mpw;
            double pzz = _mass * vz * vz * _mpw;

            for (int nx = 0; nx < SMAX; ++nx) {
                const int i = intx + nx;
                for (int ny = 0; ny < SMAX; ++ny) {
                    const int j = inty + ny;
                    for (int nz = 0; nz < SMAX; ++nz) {
                        const int k = intz + nz;
                        const auto sx = wx[nx] * wy[ny] * wz[nz];
#pragma omp atomic update
                        Pxx(i, j, k, 0) += sx * prr;
#pragma omp atomic update
                        Pyy(i, j, k, 0) += sx * ppp;
#pragma omp atomic update
                        Pzz(i, j, k, 0) += sx * pzz;
                    }
                }
            }
        }
    }
}
