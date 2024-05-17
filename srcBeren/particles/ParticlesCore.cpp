#include "containers.h"
#include "World.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "bounds.h"
#include "service.h"

void ParticlesArray::move(double dt) {
#pragma omp parallel for
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.move(dt);
        }
    }
}

void ParticlesArray::predict_velocity(const Mesh& mesh, const Domain &domain, const double dt) {

#pragma omp parallel for
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){        
            const auto coord = particle.coord;
            const auto velocity = particle.velocity;
            double3 E = get_fieldE_in_pos(mesh.fieldE, coord, domain);
            const double3 B = get_fieldB_in_pos(mesh.fieldB, coord, domain);
            const double3 En = get_fieldE_in_pos(mesh.fieldEp, coord, domain);
            E = 0.5 * (E + En);
            const auto beta = dt * charge / _mass;
            const auto alpha = 0.5*beta*mag(B);
            const auto alpha2 = alpha*alpha;
            const auto h = unit(B);
            
            const auto v12 = (1./(1.+ alpha2 )) *
                                    (velocity + alpha*cross(velocity,h) + alpha2*dot(h,velocity)*h
                                    + 0.5*beta*( E + alpha*cross(E,h) + alpha2*dot(E,h)*h) );

            particle.velocity = 2.*v12 - velocity; 
        }   
    }
}

// TO DO:: merge with void ParticlesArray::move_and_calc_current(const double
// dt, Field3d& fieldJ)
void ParticlesArray::move_and_calc_current(const double dt){

    constexpr auto SMAX = 2*SHAPE_SIZE;

    const double conx = xCellSize / (6*dt) * _mpw;
    const double cony = yCellSize / (6*dt) * _mpw;
    const double conz = zCellSize / (6*dt) * _mpw;

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        for(auto& particle : particlesData(pk)){        

                    double3 start = particle.coord;
                    
                    particle.move(dt);
                    
                    double3 end = particle.coord;

                    double xx = start.x() / xCellSize;
                    double yy = start.y() / yCellSize;
                    double zz = start.z() / zCellSize;

                    double xn = end.x() / xCellSize;
                    double yn = end.y() / yCellSize;
                    double zn = end.z() / zCellSize;
                    
                    int xk = int(xx);
                    int yk = int(yy);
                    int zk = int(zz);

                    for(int n = 0; n < SMAX; ++n){
                        for(int m = 0; m < SMAX; ++m){
                            for(int k = 0; k < SMAX; ++k){
                                jx[n][m][k] = 0.;
                                jy[n][m][k] = 0.;
                                jz[n][m][k] = 0.;
                            }
                        }
                    }
                    double arg;
                    for(int n = 0; n < SMAX; ++n){
                        arg = -xx + double(xk - GHOST_CELLS + n);
                        sx[n] = Shape2(arg);
                        arg = -yy + double(yk - GHOST_CELLS + n);
                        sy[n] = Shape2(arg);
                        arg = -zz + double(zk - GHOST_CELLS + n);
                        sz[n] = Shape2(arg);            
                        arg = -xn + double(xk - GHOST_CELLS + n);
                        sx_n[n] = Shape2(arg);
                        arg = -yn + double(yk - GHOST_CELLS + n);
                        sy_n[n] = Shape2(arg);
                        arg = -zn + double(zk - GHOST_CELLS + n);
                        sz_n[n] = Shape2(arg);
                    }

                    for(int n = 0; n < SMAX; ++n){
                        int indx = xk  + n;
                        for(int m = 0; m < SMAX; ++m){
                            int indy = yk + m;
                            for(int k = 0; k < SMAX; ++k){
                                int indz = zk + k;
                        
                                if(n == 0) jx[n][m][k] = -charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                                if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
                                
                                if(m == 0) jy[n][m][k] = -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                
                                if(k == 0) jz[n][m][k] = -charge * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-charge * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                            
                                #pragma omp atomic update
                                currentOnGrid(indx,indy,indz,0) += jx[n][m][k];
                                #pragma omp atomic update
                                currentOnGrid(indx,indy,indz,1) += jy[n][m][k];
                                #pragma omp atomic update
                                currentOnGrid(indx,indy,indz,2) += jz[n][m][k];
                            }
                        }
                    }
                }
            }

}

void ParticlesArray::move_and_calc_current(const double dt, Field3d& fieldJ){

    constexpr auto SMAX = 2*SHAPE_SIZE;

    const double conx = xCellSize / (6*dt) * _mpw;
    const double cony = yCellSize / (6*dt) * _mpw;
    const double conz = zCellSize / (6*dt) * _mpw;

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double arg;
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        for(auto& particle : particlesData(pk)){      
                    double3 start = particle.coord;
                    
                    particle.move(dt);
                    
                    double3 end = particle.coord;

                    double xx = start.x() / xCellSize;
                    double yy = start.y() / yCellSize;
                    double zz = start.z() / zCellSize;

                    double xn = end.x() / xCellSize;
                    double yn = end.y() / yCellSize;
                    double zn = end.z() / zCellSize;
                    
                    int xk = int(xx);
                    int yk = int(yy);
                    int zk = int(zz);

                    for(int n = 0; n < SMAX; ++n){
                        for(int m = 0; m < SMAX; ++m){
                            for(int k = 0; k < SMAX; ++k){
                                jx[n][m][k] = 0.;
                                jy[n][m][k] = 0.;
                                jz[n][m][k] = 0.;
                            }
                        }
                    }

                    for(int n = 0; n < SMAX; ++n){
                        arg = -xx + double(xk - GHOST_CELLS + n);
                        sx[n] = Shape2(arg);
                        arg = -yy + double(yk - GHOST_CELLS + n);
                        sy[n] = Shape2(arg);
                        arg = -zz + double(zk - GHOST_CELLS + n);
                        sz[n] = Shape2(arg);            
                        arg = -xn + double(xk - GHOST_CELLS + n);
                        sx_n[n] = Shape2(arg);
                        arg = -yn + double(yk - GHOST_CELLS + n);
                        sy_n[n] = Shape2(arg);
                        arg = -zn + double(zk - GHOST_CELLS + n);
                        sz_n[n] = Shape2(arg);
                    }

                    for(int n = 0; n < SMAX; ++n){
                        int indx = xk  + n;
                        for(int m = 0; m < SMAX; ++m){
                            int indy = yk + m;
                            for(int k = 0; k < SMAX; ++k){
                                int indz = zk + k;
                        
                                if(n == 0) jx[n][m][k] = -charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                                if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
                                
                                if(m == 0) jy[n][m][k] = -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                
                                if(k == 0) jz[n][m][k] = -charge * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-charge * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                            
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,0) += jx[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,1) += jy[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,2) += jz[n][m][k];
                            }
                        }
                    }
                }
            }
}

void ParticlesArray::calc_Esirkepov_current(const double dt, Field3d& fieldJ) const{
    constexpr auto SMAX = 2*SHAPE_SIZE;
        const double conx = xCellSize / (6*dt) * _mpw;
        const double cony = yCellSize / (6*dt) * _mpw;
        const double conz = zCellSize / (6*dt) * _mpw;
#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double arg;
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        double q = charge;

        for(auto& particle : particlesData(pk)){      
                    const auto start = particle.initCoord;                    
                    const auto end = particle.coord;

                    const auto xx = start.x() / xCellSize;
                    const auto yy = start.y() / yCellSize;
                    const auto zz = start.z() / zCellSize;

                    const auto xn = end.x() / xCellSize;
                    const auto yn = end.y() / yCellSize;
                    const auto zn = end.z() / zCellSize;
                    
                    const auto xk = int(xx);
                    const auto yk = int(yy);
                    const auto zk = int(zz);

                    for(int n = 0; n < SMAX; ++n){
                        for(int m = 0; m < SMAX; ++m){
                            for(int k = 0; k < SMAX; ++k){
                                jx[n][m][k] = 0.;
                                jy[n][m][k] = 0.;
                                jz[n][m][k] = 0.;
                            }
                        }
                    }

                    for(int n = 0; n < SMAX; ++n){
                        arg = -xx + double(xk - GHOST_CELLS + n);
                        sx[n] = Shape(arg);
                        arg = -yy + double(yk - GHOST_CELLS + n);
                        sy[n] = Shape(arg);
                        arg = -zz + double(zk - GHOST_CELLS + n);
                        sz[n] = Shape(arg);            
                        arg = -xn + double(xk - GHOST_CELLS + n);
                        sx_n[n] = Shape(arg);
                        arg = -yn + double(yk - GHOST_CELLS + n);
                        sy_n[n] = Shape(arg);
                        arg = -zn + double(zk - GHOST_CELLS + n);
                        sz_n[n] = Shape(arg);
                    }

                    for(int n = 0; n < SMAX; ++n){
                        const auto indx = xk  + n;
                        for(int m = 0; m < SMAX; ++m){
                            const auto indy = yk + m;
                            for(int k = 0; k < SMAX; ++k){
                                const auto indz = zk + k;
                        
                                if(n == 0) jx[n][m][k] = -q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                                if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
                                
                                if(m == 0) jy[n][m][k] = -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                
                                if(k == 0) jz[n][m][k] = -q * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-q * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                #pragma omp atomic update                            
                                fieldJ(indx,indy,indz,0) += jx[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,1) += jy[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,2) += jz[n][m][k];
                            }
                        }
                    }
                }
            }
}

void ParticlesArray::correctv(Mesh& mesh, const Domain& domain, const double dt){

    std::array<double,20> ldistr;
    for (auto& val: ldistr){
        val =0.0;
    }

    double jp_cell = 0;
#pragma omp parallel for reduction(+:jp_cell)
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto initVelocity = particle.initVelocity;
                    const auto velocity = particle.velocity;
 
                    double3 end = particle.coord;
                    double3 coord = end - 0.5*dt*velocity;
                    double3 Ep = get_fieldE_in_pos(mesh.fieldEp, coord, domain);
                    double3 E = get_fieldE_in_pos(mesh.fieldE, coord, domain);

                    double3 v12 = 0.5*(velocity + initVelocity); 

                    jp_cell += 0.5*_mpw*charge*dot(v12,(Ep+E));
                }
            }

    mesh.fieldJp_full.data() = mesh.fieldJp.data() + mesh.Lmat2*(mesh.fieldE.data() + mesh.fieldEp.data())/dt;

    const double energyJeEn = mesh.calc_JE(mesh.fieldEn,currentOnGrid);
    const double energyJeE = mesh.calc_JE(mesh.fieldE,currentOnGrid);
    const double energyJpEp = mesh.calc_JE(mesh.fieldEp,mesh.fieldJp_full);
    const double energyJpE = mesh.calc_JE(mesh.fieldE,mesh.fieldJp_full);
    const double energyK = get_kinetic_energy();
    const double lambda = sqrt(1 + dt*( 0.5*(energyJeEn + energyJeE) - jp_cell) / energyK);

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto velocity = particle.velocity;

                    particle.velocity = lambda*velocity;
                }
            }

    const double energyK2 = get_kinetic_energy();
    std::cout << "lambda "<< lambda  << " " << lambda*lambda << " "<< energyK2-energyK << " " << 0.5*dt*(energyJeEn + energyJeE - energyJpEp - energyJpE) << "\n";

}


void ParticlesArray::correctv_component(Mesh& mesh, const Domain &domain, const double dt){
    double jp_cellx = 0;
    double jp_celly = 0;
    double jp_cellz = 0;
#pragma omp parallel for reduction(+ : jp_cellx) reduction(+ : jp_celly) reduction(+:jp_cellz)
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto initVelocity = particle.initVelocity;
                    const auto velocity = particle.velocity;
 
                    double3 end = particle.coord;
                    double3 coord = end - 0.5*dt*velocity;
                    
                    double3 Ep = get_fieldE_in_pos(mesh.fieldEp,coord, domain); 
                    double3 E = get_fieldE_in_pos(mesh.fieldE,coord, domain);
                    E+= Ep;

                    double3 v12 = 0.5*(velocity + initVelocity);
                    double3 vE = double3(v12.x() * E.x(), v12.y() * E.y(),
                                         v12.z() * E.z() );

                    jp_cellx += (0.5 * _mpw * charge) * vE.x();
                    jp_celly += (0.5 * _mpw * charge) * vE.y();
                    jp_cellz += (0.5 * _mpw * charge) * vE.z();
        }
            }


    const double3 energyJeEn = mesh.calc_JE_component(mesh.fieldEn,currentOnGrid);
    const double3 energyJeE =
        mesh.calc_JE_component(mesh.fieldE, currentOnGrid);
    const double3 energyK = get_kinetic_energy_component();
    double3 lambda;
    lambda.x() = 
            sqrt(1 + dt * (0.5 * (energyJeEn.x() + energyJeE.x()) - jp_cellx) /
                     energyK.x());
    lambda.y() =
          sqrt(1 + dt * (0.5 * (energyJeEn.y() + energyJeE.y()) - jp_celly) /
                     energyK.y());
    lambda.z() =
         sqrt(1 + dt * (0.5 * (energyJeEn.z() + energyJeE.z()) - jp_cellz) /
                     energyK.z());
    // double lambda2 =
    //     sqrt(1 + Dt *
    //                  (0.5 * (energyJeEn.x() + energyJeE.x() + energyJeEn.y() +
    //                          energyJeE.y() + energyJeEn.z() + energyJeE.z()) -
    //                   jp_cellx - jp_celly - jp_cellz) /
    //                  (energyK.x() + energyK.y() + energyK.z()));

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto velocity = particle.velocity;

                    particle.velocity.x() = lambda.x() * velocity.x();
                    particle.velocity.y() = lambda.y() * velocity.y();
                    particle.velocity.z() = lambda.z() * velocity.z();
                    //particle.velocity = lambda2 * velocity;
        }
            }

    //const double energyK2 = get_kinetic_energy();
    std::cout << "lambda "<< lambda << "\n";

}

void ParticlesArray::predict_current(const Field3d& fieldB, Field3d& fieldJ,
                                     const Domain& domain, const double dt) {
    constexpr auto SMAX = SHAPE_SIZE; 
#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        int i,j,k;
        int i05,j05,k05;
        double qp = charge;

        alignas(64) double wx[SMAX], wy[SMAX], wz[SMAX];
        alignas(64) double wx05[SMAX], wy05[SMAX], wz05[SMAX];

        for(auto& particle : particlesData(pk)){     
                    double3 coord = particle.coord;
                    double3 velocity = particle.velocity;

                    double x = coord.x() / xCellSize + GHOST_CELLS;
                    double y = coord.y() / yCellSize + GHOST_CELLS;
                    double z = coord.z() / zCellSize + GHOST_CELLS;
                    double x05 = x - 0.5;
                    double y05 = y - 0.5;
                    double z05 = z - 0.5;

                    const auto ix = int(x);
                    const auto iy = int(y);
                    const auto iz = int(z);
                    const auto ix05 = int(x05);
                    const auto iy05 = int(y05);
                    const auto iz05 = int(z05);
                    
                    wx[1] = (x - ix);
                    wx[0] = 1 - wx[1];
                    wy[1] = (y - iy);
                    wy[0] = 1 - wy[1];
                    wz[1] = (z - iz);
                    wz[0] = 1 - wz[1];

                    wx05[1] = (x05 - ix05);
                    wx05[0] = 1 - wx05[1];
                    wy05[1] = (y05 - iy05);
                    wy05[0] = 1 - wy05[1];
                    wz05[1] = (z05 - iz05);
                    wz05[0] = 1 - wz05[1];        
                    
                    double3 B = get_fieldB_in_pos(fieldB,coord, domain); 

                    double beta = dt * qp / _mass;
                    double alpha = 0.5*beta*mag(B);
                    double alpha2 = alpha*alpha;
                    double3 h = unit(B);

                    double3 current = qp * _mpw / (1.+alpha2) * (velocity + alpha*cross(velocity,h) + alpha2*dot(h,velocity)*h );

                    for(int nx = 0; nx < SMAX; ++nx){
                        i = ix + nx;
                        i05 = ix05 + nx;
                        for(int ny = 0; ny < SMAX; ++ny){
                            j = iy  + ny;
                            j05 = iy05  + ny;
                            for(int nz = 0; nz < SMAX; ++nz){
                                k = iz  + nz;
                                k05 = iz05  + nz;
                                double sx = wx05[nx] * wy[ny] * wz[nz];
                                double sy = wx[nx] * wy05[ny] * wz[nz];
                                double sz = wx[nx] * wy[ny] * wz05[nz];
                                #pragma omp atomic update
                                fieldJ(i05,j,k,0) += sx*current.x();
                                #pragma omp atomic update
                                fieldJ(i,j05,k,1) += sy*current.y();
                                #pragma omp atomic update
                                fieldJ(i,j,k05,2) += sz*current.z();                   
                            }
                        }
                    }

                }
            }
}


// Very slow function. Fill Lmatrix by each particles
void ParticlesArray::get_L(Mesh& mesh, const Domain &domain, const double dt) {
#pragma omp parallel
{
    for (int xStep = 0; xStep < 4; xStep++) {
        for (int yStep = 0; yStep < 4; yStep++) {
            #pragma omp for collapse(2)
            for (int ix = xStep; ix < particlesData.size().x(); ix += 4) {
                for (int iy = yStep; iy < particlesData.size().y(); iy += 4) {
                    for (int iz = 0; iz < particlesData.size().z(); ++iz) {
                        for (auto& particle : particlesData(ix,iy,iz)) {
                            const auto coord = particle.coord;
                            mesh.update_Lmat(coord, domain, charge, _mass, _mpw, dt);
                        }
                    }
                }
            }
#pragma omp barrier
        }
    }
}
}

