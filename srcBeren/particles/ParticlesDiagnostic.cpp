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

// void ParticlesArray::glue_density_bound()
// { 
    
//     int overlap = 3;
//     int3 size = densityOnGrid.size();
//     // TODO :: add ccheck for periodic 
//     //if(isPeriodicX){
//         for (auto i = 0; i < overlap; ++i){
//           for (auto l = 0; l < size.y(); ++l){
//             for (auto k = 0; k < size.z(); ++k){
//                 densityOnGrid(i,l,k) += densityOnGrid(i + size.x() - overlap,l,k);
//                 densityOnGrid(i + size.x() - overlap,l,k) = densityOnGrid(i,l,k);
//             }
//           }
//         }
//     //}
//     //if(isPeriodicY){
//         for (auto i = 0; i < size.x(); ++i){
//           for (auto l = 0; l < overlap; ++l){
//             for (auto k = 0; k < size.z(); ++k){
//                 densityOnGrid(i,l,k) += densityOnGrid(i,l + size.y() - overlap,k);
//                 densityOnGrid(i,l+ size.y() - overlap,k) = densityOnGrid(i,l ,k);
//             }
//           }
//         }
//     //}
//     //if(isPeriodicZ){
//         for (auto i = 0; i < size.x(); ++i){
//           for (auto l = 0; l < size.y(); ++l){
//             for (auto k = 0; k < overlap; ++k){
//                 densityOnGrid(i,l,k) += densityOnGrid(i ,l,k+ size.z() - overlap);
//                 densityOnGrid(i ,l,k+ size.z() - overlap) = densityOnGrid(i,l,k);
//             }
//           }
//         }
//     //}
// }

void ParticlesArray::density_on_grid_update(){

    auto ShapeN = Shape;

     constexpr auto SMAX = 2 * SHAPE_SIZE;

    densityOnGrid.set_zero();

#pragma omp parallel for
    for ( auto j = 0; j < size(); ++j){
        double arg;
    alignas(64) double sx[SMAX];
    alignas(64) double sy[SMAX];
    alignas(64) double sz[SMAX];

        for(const auto& particle : particlesData(j)){
            int xk = int(particle.coord.x() / xCellSize);
            int yk = int(particle.coord.y() / yCellSize);
            int zk = int(particle.coord.z() / zCellSize);

            for(int n = 0; n < SMAX; ++n){
                arg = -particle.coord.x() / xCellSize + double(xk - GHOST_CELLS + n);
                sx[n] = ShapeN(arg);
                arg = -particle.coord.y() / yCellSize +
                      double(yk - GHOST_CELLS + n);
                sy[n] = ShapeN(arg);
                arg = -particle.coord.z() / zCellSize +
                      double(zk - GHOST_CELLS + n);
                sz[n] = ShapeN(arg);
            }

            for(int n = 0; n < SMAX; ++n){
                int indx = xk + n;
                for(int m = 0; m < SMAX; ++m){
                    int indy = yk + m;
                    for(int k = 0; k < SMAX; ++k){
                    int indz = zk + k;
#pragma omp atomic update
                    densityOnGrid(indx,indy,indz,0) += _mpw*charge * sx[n] * sy[m] * sz[k];
                    } 
                }
            }
        }
    }
    make_periodic_border_with_add(densityOnGrid, bounds);
}

void ParticlesArray::phase_on_grid_update(const Domain& domain){ 
    // const double3 cellSize = domain.cell_size();
    // int pk, xk;
    // double x, px;
    // bool blounder;
    // Particle particle;
    // double x_min = 0.;
    // double x_max = cellSize.x()*domain.num_cells(Dim::X);
    // double pdx = (x_max - x_min) / 100; // TO DO: delete magic numbers
    // double pdp =
    //     (phasePXmax - phasePXmin) / 100;   // TO DO: delete magic numbers

    // phaseOnGrid.set_zero();

    // for (auto k = 0 ; k < size(); k++){
    //     for(const auto& particle : particlesData(k)){
    //         x = particle.coord.x();
    //         px = particle.velocity.x();

    //         xk = int((x-x_min) / pdx); 
    //         pk = int((px - phasePXmin) / pdp);

    //         blounder = (xk < 0) || (xk >= 100) || (pk < 0) ||
    //                    (pk >= 100);   // TO DO: delete magic numbers
    //         if(!blounder){
    //             phaseOnGrid(xk, pk) +=
    //                 (mpw() / (cellSize.x() * cellSize.y() * cellSize.z() * pdx * pdp));
    //         }
    //     }
    // }
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
    make_periodic_border_with_add(Pxx, bounds);
    make_periodic_border_with_add(Pyy, bounds);
    make_periodic_border_with_add(Pzz, bounds);
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
    make_periodic_border_with_add(Pxx, bounds);
    make_periodic_border_with_add(Pyy, bounds);
    make_periodic_border_with_add(Pzz, bounds);
}
