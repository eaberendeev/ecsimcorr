#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "solverSLE.h"
#include "util.h"

void Mesh::init(const Domain &domain, const ParametersMap &parameters){
    Lmat.resize(domain.total_size() * 3, domain.total_size() * 3);
        Lmat2.resize(domain.total_size() * 3,
                     domain.total_size() * 3);
        Mmat.resize(domain.total_size() * 3,
                    domain.total_size() * 3);
        Imat.resize(domain.total_size() * 3,
                    domain.total_size() * 3);
        curlE.resize(domain.total_size() * 3,
                     domain.total_size() * 3);
        curlB.resize(domain.total_size() * 3,
                     domain.total_size() * 3);

        LmatX.resize(domain.total_size() * 3);
        chargeDensityOld.resize(domain.size(), 1);
        chargeDensity.resize(domain.size(), 1);
        divE.resize(domain.total_size(), domain.total_size() * 3);

        xCellSize = domain.cell_size().x();
        yCellSize = domain.cell_size().y();
        zCellSize = domain.cell_size().z();
        xSize = domain.size().x();
        ySize = domain.size().y();
        zSize = domain.size().z();

        stencil_Imat(domain);
        stencil_curlE(domain);
        stencil_curlB(domain);
        stencil_divE(domain);
        bounds.setBounds(domain.lower_bounds(), domain.upper_bounds());

        double dt = parameters.get_double("Dt");
        Mmat = -0.25 * dt * dt * curlB * curlE;
        for (auto& rowMap : LmatX) {
            rowMap.reserve(LMAT_MAX_ELEMENTS_PER_ROW);
        }
}

void Mesh::print_operator(const Operator &oper){
     for (int k=0; k < oper.outerSize(); ++k){
       for (Eigen::SparseMatrix<double,MAJOR>::InnerIterator it(oper,k); it; ++it){
         std::cout  << pos_vind(it.row(),0) << " " 
         << pos_vind(it.row(),1) << " "
         << pos_vind(it.row(),2) << " " 
         << pos_vind(it.row(),3) << " " 
         <<  pos_vind(it.col(),0) << " " 
         <<  pos_vind(it.col(),1) << " " 
         <<  pos_vind(it.col(),2) << " " 
         <<  pos_vind(it.col(),3) << " " << it.value()  << "\n";
       }
    }
}

void Mesh::set_uniform_field(Field3d& field, double bx, double by, double bz) {
    auto size = field.size();
    for (auto i = 0; i < size.x(); ++i) {
        for (auto j = 0; j < size.y(); ++j) {
            for (auto k = 0; k < size.z(); ++k) {
                field(i, j, k, Dim::X) = bx;
                field(i, j, k, Dim::Y) = by;
                field(i, j, k, Dim::Z) = bz;
            }
        }
    }
}

void Mesh::prepare()
{
//   fieldE = fieldEn;
//   fieldB0 = fieldB;
//   fieldB = fieldBn;

  //Lmat2.setZero();

#pragma omp parallel for
  for ( size_t i = 0; i < LmatX.size(); i++){
      for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
        it->second = 0.;
    }
  }
}

double Mesh::calc_energy_field(const Field3d& field) const{
  double potE = 0;
  int i_max = field.size().x();
  int j_max = field.size().y();
  int k_max = field.size().z();
  constexpr int OVERLAP_SIZE = 3;
  if (bounds.isPeriodic(X)) {
      i_max -= OVERLAP_SIZE;
  }
  if (bounds.isPeriodic(Y)) {
      j_max -= OVERLAP_SIZE;
  }
  if (bounds.isPeriodic(Z)) {
      k_max -= OVERLAP_SIZE;
  }

  for(auto i = 0; i < i_max; ++i){
    for(auto j = 0; j < j_max; ++j){
      for(auto k = 0; k < k_max; ++k){
        double3 v = double3(field(i,j,k,0),field(i,j,k,1),field(i,j,k,2));
        potE += dot(v,v );      
      }
    }
  }
  potE *= (0.5);
  return potE;
}

double calc_JE(const Field3d& fieldE,const Field3d& fieldJ, const Bounds& bounds) {
  double potE = 0;
  int i_max = fieldE.size().x();
  int j_max = fieldE.size().y();
  int k_max = fieldE.size().z();

  constexpr int OVERLAP_SIZE = 3;
  if (bounds.isPeriodic(X)) {
      i_max -= OVERLAP_SIZE;
  }
  if (bounds.isPeriodic(Y)) {
      j_max -= OVERLAP_SIZE;
  }
  if (bounds.isPeriodic(Z)) {
      k_max -= OVERLAP_SIZE;
  }
  for(auto i = 0; i < i_max; ++i){
    for(auto j = 0; j < j_max; ++j){
      for(auto k = 0; k < k_max; ++k){
        double3 E = double3(fieldE(i,j,k,0),fieldE(i,j,k,1),fieldE(i,j,k,2));
        double3 J = double3(fieldJ(i,j,k,0),fieldJ(i,j,k,1),fieldJ(i,j,k,2));
        potE += dot(J,E);      
      }
    }
  }
  return potE;
}

double3 calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ, const Bounds& bounds) {
    double3 potE = double3(0,0,0);
    int i_max = fieldE.size().x();
    int j_max = fieldE.size().y();
    int k_max = fieldE.size().z();

    constexpr int OVERLAP_SIZE = 3;
    if (bounds.isPeriodic(X)) {
        i_max -= OVERLAP_SIZE;
    }
    if (bounds.isPeriodic(Y)) {
        j_max -= OVERLAP_SIZE;
    }
    if (bounds.isPeriodic(Z)) {
        k_max -= OVERLAP_SIZE;
    }
    for (auto i = 0; i < i_max; ++i) {
        for (auto j = 0; j < j_max; ++j) {
            for (auto k = 0; k < k_max; ++k) {
                double3 E = double3(fieldE(i, j, k, 0), fieldE(i, j, k, 1),
                                    fieldE(i, j, k, 2));
                double3 J = double3(fieldJ(i, j, k, 0), fieldJ(i, j, k, 1),
                                    fieldJ(i, j, k, 2));
                potE += double3(J.x() * E.x(), J.y() * E.y(), J.z() * E.z() );
            }
        }
    }
    return potE;
}

// Solve Ax=b for find fieldE
void Mesh::correctE(Field3d& En, const Field3d& E, const Field3d& B,
                    const Field3d& J, const double dt) {
    // get x
	Field rhs = E.data() - dt*J.data() + dt*curlB*B.data() + Mmat*E.data();
    Operator A = Imat - Mmat;
    // solve Ax=b, fieldEn - output
    solve_SLE(A, rhs, En.data(), E.data());
    
    std::cout<< "Solver2 error = "<< (A*En.data() - rhs).norm() << "\n";
}

// Solve Ax=b for find fieldE.
// (E_{n+1} - E_n) / dt = -J_{n+1/2} + rot(B_{n+1/2}) B_{n+1/2} =
// (B_n + B_{n+1})/2 (B_{n+1}- B_n) / dt = - rot(E_{n+1/2})
// M = -0.25 * dt * dt * rot_opB * rot_opE;
// E_{n+1} = dt*E_n + M*(E_{n+1}+E_n) - dt*J_{n+1/2}+ rot(B_n)

// E_{n+1} - fieldEnew (out)
// E_n - fieldE (in)
// B_n - fieldB (in)
// J_{n+1/2} - fieldJ (in)
void Mesh::impicit_find_fieldE(Field3d& Enew, const Field3d& E, const Field3d& B,
                    const Field3d& J, const double dt) {

    Field rhs = E.data() - dt * J.data() +
                dt * curlB * B.data() + Mmat * E.data();
    Operator A = Imat - Mmat;

    solve_SLE(A, rhs, Enew.data(), E.data());

    std::cout << "Solver impicit_find_fieldE error = "
              << (A * Enew.data() - rhs).norm() << "\n";
}

double Mesh::calculate_residual(const Field3d& Enew, const Field3d& E,
                               const Field3d& B, const Field3d& J,
                               const double dt) {
    Field rhs =
        E.data() - dt * J.data() + dt * curlB * B.data() + Mmat * E.data();
    Operator A = Imat - Mmat;

    return (A * Enew.data() - rhs).norm();
}

void Mesh::predictE(Field3d& Ep, const Field3d& E, const Field3d& B,
                    const Field3d& J, double dt) {
    Operator Lmat2 = Mmat - Lmat;

    Field rhs =
        E.data() - dt * J.data() + dt * curlB * B.data() + Lmat2 * E.data();
    Operator A = Imat - Lmat2;
    solve_SLE(A, rhs, Ep.data(), E.data());

    std::cout << "Solver1 error = " << (A * Ep.data() - rhs).norm()
              << "\n";
}

void Mesh::fdtd_explicit(Field3d& E, Field3d& B, const Field3d& J,
                         const double dt) {
    E.data() +=
        0.5 * dt * curlB * B.data() - 0.5 * dt * J.data();
    B.data() -= 0.5 * dt * curlE * E.data();
}

void Mesh::computeB(const Field3d& fieldE, const Field3d& fieldEn,
                    Field3d& fieldB, double dt) {
    fieldB.data() -= 0.5*dt*curlE*(fieldE.data() + fieldEn.data() );
}

void Mesh::compute_fieldB(Field3d& Bn, const Field3d& B, const Field3d& E,
                          const Field3d& En, double dt){
    Bn.data() = B.data() - 0.5 * dt * curlE * (E.data() + En.data());
}
// std::tuple<int, int, int> calculateCellIndices(const double3& coord,
//                                                const Domain& domain) {
//     int cellIndexX = coord.x() / domain.cellSize.x + domain.GHOST_CELLS;
//     int cellIndexY = coord.y() / domain.Dy + domain.GHOST_CELLS;
//     int cellIndexZ = coord.z() / domain.Dz + domain.GHOST_CELLS;

//     return {cellIndexX, cellIndexY, cellIndexZ};
// }

// std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
// calculateShapeCoefficients(const double3& coord, int cellIndexX, int cellIndexY,
//                            int cellIndexZ) {
//     //...
//     return {sx, sy, sz};
// }

// void update_Lmat(const double3& coord, const Domain& domain, double charge,
//                  double mass, double mpw, double dt) {
//     const int SMAX = SHAPE_SIZE;

//     // Вынесли расчет индексов ячеек в отдельную функцию
//     auto [cellIndexX, cellIndexY, cellIndexZ] =
//         calculateCellIndices(coord, domain);

//     // Вынесли расчет коэффициентов формы в отдельную функцию
//     auto [sx, sy, sz] =
//         calculateShapeCoefficients(coord, cellIndexX, cellIndexY, cellIndexZ);

//     auto B = calculateMagneticField(cellIndexX, cellIndexY, cellIndexZ, sx, sy,
//                                     sz, domain);

//     // Проверка границ цикла
//     for (int i = 0; i < SMAX && i < sx.size(); ++i) {
//         for (int j = 0; j < SMAX && j < sy.size(); ++j) {
//             for (int k = 0; k < SMAX && k < sz.size(); ++k) {
//                 auto value = calculateLmatValue(i, j, k, sx, sy, sz, B, charge,
//                                                 mass, mpw, dt);

//                 if (std::fabs(value) > 1e-16) {
//                     auto indx1 = vind(/*...*/);
//                     auto indx2 = vind(/*...*/);
//                     LmatX[indx1][indx2] += value;
//                 }
//             }
//         }
//     }
// }


void Mesh::update_Lmat(const double3& coord, const Domain &domain, double charge, double mass,
                       double mpw, const Field3d& fieldB, const double dt) {
    const int SMAX = SHAPE_SIZE;
    double3 B;
    double wx,wy,wz,wx1,wy1,wz1;
    double value;
    int cellLocX,cellLocY,cellLocZ,cellLocX05,cellLocY05,cellLocZ05;
    double coordLocX,coordLocY,coordLocZ;
    double coordLocX05,coordLocY05,coordLocZ05;
    int i,j,k,i1,j1,k1;
    int indx1,indy1,indz1;
    int indx2,indy2,indz2;
    int indx,indy,indz;
    int indx05,indy05,indz05;
    alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
    alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];

        B = 0.;

        coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
        coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
        coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
        coordLocX05 = coordLocX - 0.5;
        coordLocY05 = coordLocY - 0.5;
        coordLocZ05 = coordLocZ - 0.5;

        cellLocX = int(coordLocX);
        cellLocY = int(coordLocY);
        cellLocZ = int(coordLocZ);
        cellLocX05 = int(coordLocX05);
        cellLocY05 = int(coordLocY05);
        cellLocZ05 = int(coordLocZ05);
        
        sx[1] = (coordLocX - cellLocX);
        sx[0] = 1 - sx[1];
        sy[1] = (coordLocY - cellLocY);
        sy[0] = 1 - sy[1];
        sz[1] = (coordLocZ - cellLocZ);
        sz[0] = 1 - sz[1];

        sx05[1] = (coordLocX05 - cellLocX05);
        sx05[0] = 1 - sx05[1];
        sy05[1] = (coordLocY05 - cellLocY05);
        sy05[0] = 1 - sy05[1];
        sz05[1] = (coordLocZ05 - cellLocZ05);
        sz05[0] = 1 - sz05[1];        

        for(i = 0; i < SMAX; ++i){
            indx = cellLocX + i;
            indx05 = cellLocX05 + i;
            for(j = 0; j < SMAX; ++j){
                indy = cellLocY  + j;
                indy05 = cellLocY05  + j;
                for(k = 0; k < SMAX; ++k){
                    indz = cellLocZ  + k;
                    indz05 = cellLocZ05  + k;
                    wx = sx[i] * sy05[j] * sz05[k];
                    wy = sx05[i] * sy[j] * sz05[k];
                    wz = sx05[i] * sy05[j] * sz[k];
                    B.x() += (wx * fieldB(indx,indy05,indz05,0) );
                    B.y() += (wy * fieldB(indx05,indy,indz05,1) );
                    B.z() += (wz * fieldB(indx05,indy05,indz,2) );
                }
            }
        }
        const double3 h = unit(B);

        const double alpha = 0.5*dt*charge*mag(B) / mass;
        const double Ap = 0.25*dt*dt*mpw*charge*charge / mass / (1+alpha*alpha); 
        
        for(i = 0; i < SMAX; ++i){
            indx = cellLocX + i;
            indx05 = cellLocX05 + i;
            for(j = 0; j < SMAX; ++j){
                indy = cellLocY  + j;
                indy05 = cellLocY05  + j;
                for(k = 0; k < SMAX; ++k){
                    indz = cellLocZ  + k;
                    indz05 = cellLocZ05  + k;
                    wx = sx05[i] * sy[j] * sz[k];
                    wy = sx[i] * sy05[j] * sz[k];
                    wz = sx[i] * sy[j] * sz05[k];
                    
                    for(i1 = 0; i1 < SMAX; ++i1){
                        for(j1 = 0; j1 < SMAX; ++j1){
                            for(k1 = 0; k1 < SMAX; ++k1){
                                wx1 = sx05[i1] * sy[j1] * sz[k1];
                                wy1 = sx[i1] * sy05[j1] * sz[k1];
                                wz1 = sx[i1] * sy[j1] * sz05[k1];
                                // xx
                                value = wx*wx1*Ap*(1.+alpha*alpha*h.x()*h.x() );
                                if(fabs(value) > 1.e-16){
                                    indx1 = vind(cellLocX05 + i,cellLocY + j,cellLocZ+k,0);
                                    indx2 = vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0);
                                    LmatX[indx1][indx2] += value;
                                }
                                // xy
                                value = wx*wy1*Ap*(alpha*h.z() + alpha*alpha*h.x()*h.y() );
                                if(fabs(value) > 1.e-16){
                                    indx1 = vind(cellLocX05 + i,cellLocY + j,cellLocZ+k,0);
                                    indy2 = vind(cellLocX + i1,cellLocY05 + j1,cellLocZ+k1,1);
                                    LmatX[indx1][indy2] += value;
                                }
                                 // xz
                                  value = wx*wz1*Ap*alpha*(-h.y() + alpha*h.x() *h.z() );
                                  if(fabs(value) > 1.e-16){
                                    indx1 = vind(cellLocX05 + i,cellLocY + j,cellLocZ+k,0);
                                    indz2 = vind(cellLocX + i1,cellLocY + j1,cellLocZ05+k1,2);
                                    LmatX[indx1][indz2] += value;                                  
                                  }
                                // yx
                                value = wy*wx1*Ap*alpha*(-h.z() + alpha*h.x() *h.y() );
                                 if(fabs(value) > 1.e-16){
                                    indy1 = vind(cellLocX + i,cellLocY05 + j,cellLocZ+k,1);
                                    indx2 = vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0);
                                    LmatX[indy1][indx2] += value;
                                  } 
                                // yy 
                                value = wy*wy1*Ap*(1.+alpha*alpha*h.y()*h.y() );
                                if(fabs(value) > 1.e-16){
                                    indy1 = vind(cellLocX + i,cellLocY05 + j,cellLocZ+k,1);
                                    indy2 = vind(cellLocX + i1,cellLocY05 + j1,cellLocZ+k1,1);
                                    LmatX[indy1][indy2] += value;
                                }
                                // yz
                                value = wy*wz1*Ap*alpha*(h.x() + alpha*h.y() *h.z() );
                                if(fabs(value) > 1.e-16){
                                    indy1 = vind(cellLocX + i,cellLocY05 + j,cellLocZ+k,1);
                                    indz2 = vind(cellLocX + i1,cellLocY + j1,cellLocZ05+k1,2);
                                    LmatX[indy1][indz2] += value;
                                }
                                // zx
                                value = wz*wx1*Ap*alpha*(h.y() + alpha*h.x() *h.z() );
                                if(fabs(value) > 1.e-16){
                                    indz1 = vind(cellLocX + i,cellLocY + j,cellLocZ05+k,2);
                                    indx2 = vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0);
                                    LmatX[indz1][indx2] += value;
                                }
                                value = wz*wy1*Ap*alpha*(-h.x() + alpha*h.y() *h.z() );
                                if(fabs(value) > 1.e-16){
                                    indz1 = vind(cellLocX + i,cellLocY + j,cellLocZ05+k,2);
                                    indy2 = vind(cellLocX + i1,cellLocY05 + j1,cellLocZ+k1,1);
                                    LmatX[indz1][indy2] += value;
                                }
                                value = wz*wz1*Ap*(1.+alpha*alpha*h.z()*h.z() );
                                if(fabs(value) > 1.e-16){
                                    indz1 = vind(cellLocX + i,cellLocY + j,cellLocZ05+k,2);
                                    indz2 = vind(cellLocX + i1,cellLocY + j1,cellLocZ05+k1,2);
                                    LmatX[indz1][indz2] += value;
                                }
                            } // G'z
                        } // G'y
                    } // G'x
                } // Gz
            } // Gy
        } // Gx
}

void Mesh::update_Lmat2(const double3& coord, const Domain& domain,
                       double charge, double mass, double mpw,
                       const Field3d& fieldB, const double dt) {
    const int SMAX = SHAPE_SIZE;
    double3 B;
    double wx, wy, wz, wx1, wy1, wz1;
    double value;
    int cellLocX, cellLocY, cellLocZ, cellLocX05, cellLocY05, cellLocZ05;
    double coordLocX, coordLocY, coordLocZ;
    double coordLocX05, coordLocY05, coordLocZ05;
    int i, j, k, i1, j1, k1;
    int indx1, indy1, indz1;
    int indx2, indy2, indz2;
    int indx, indy, indz;
    int indx05, indy05, indz05;
    alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
    alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];

    B = 0.;

    coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    coordLocX05 = coordLocX - 0.5;
    coordLocY05 = coordLocY - 0.5;
    coordLocZ05 = coordLocZ - 0.5;

    cellLocX = int(coordLocX);
    cellLocY = int(coordLocY);
    cellLocZ = int(coordLocZ);
    cellLocX05 = int(coordLocX05);
    cellLocY05 = int(coordLocY05);
    cellLocZ05 = int(coordLocZ05);

    sx[1] = (coordLocX - cellLocX);
    sx[0] = 1 - sx[1];
    sy[1] = (coordLocY - cellLocY);
    sy[0] = 1 - sy[1];
    sz[1] = (coordLocZ - cellLocZ);
    sz[0] = 1 - sz[1];

    sx05[1] = (coordLocX05 - cellLocX05);
    sx05[0] = 1 - sx05[1];
    sy05[1] = (coordLocY05 - cellLocY05);
    sy05[0] = 1 - sy05[1];
    sz05[1] = (coordLocZ05 - cellLocZ05);
    sz05[0] = 1 - sz05[1];

    for (i = 0; i < SMAX; ++i) {
        indx = cellLocX + i;
        indx05 = cellLocX05 + i;
        for (j = 0; j < SMAX; ++j) {
            indy = cellLocY + j;
            indy05 = cellLocY05 + j;
            for (k = 0; k < SMAX; ++k) {
                indz = cellLocZ + k;
                indz05 = cellLocZ05 + k;
                wx = sx[i] * sy05[j] * sz05[k];
                wy = sx05[i] * sy[j] * sz05[k];
                wz = sx05[i] * sy05[j] * sz[k];
                B.x() += (wx * fieldB(indx, indy05, indz05, 0));
                B.y() += (wy * fieldB(indx05, indy, indz05, 1));
                B.z() += (wz * fieldB(indx05, indy05, indz, 2));
            }
        }
    }
    const double3 h = unit(B);

    const double alpha = 0.5 * dt * charge * mag(B) / mass;
    const double Ap =
        0.25 * dt * dt * mpw * charge * charge / mass / (1 + alpha * alpha);

    for (i = 0; i < SMAX; ++i) {
        indx = cellLocX + i;
        indx05 = cellLocX05 + i;
        for (j = 0; j < SMAX; ++j) {
            indy = cellLocY + j;
            indy05 = cellLocY05 + j;
            for (k = 0; k < SMAX; ++k) {
                indz = cellLocZ + k;
                indz05 = cellLocZ05 + k;
                wx = sx05[i] * sy[j] * sz[k];
                wy = sx[i] * sy05[j] * sz[k];
                wz = sx[i] * sy[j] * sz05[k];

                for (i1 = 0; i1 < SMAX; ++i1) {
                    for (j1 = 0; j1 < SMAX; ++j1) {
                        for (k1 = 0; k1 < SMAX; ++k1) {
                            wx1 = sx05[i1] * sy[j1] * sz[k1];
                            wy1 = sx[i1] * sy05[j1] * sz[k1];
                            wz1 = sx[i1] * sy[j1] * sz05[k1];
                            // xx
                            value = wx * wx1 * Ap *
                                    (1. + alpha * alpha * h.x() * h.x());

                                indx1 = vind(cellLocX05 + i, cellLocY + j,
                                             cellLocZ + k, 0);
                                // indx2 = vind(cellLocX05 + i1, cellLocY + j1,
                                //              cellLocZ + k1, 0);
                                indx2 = get_col_index_Lx(i1,
                                                         j1,
                                                         k1, X);
                                LmatX2[indx1][indx2] += value;

                            // xy
                            value =
                                wx * wy1 * Ap *
                                (alpha * h.z() + alpha * alpha * h.x() * h.y());

                                indx1 = vind(cellLocX05 + i, cellLocY + j,
                                             cellLocZ + k, 0);
                                // indy2 = vind(cellLocX + i1, cellLocY05 + j1,
                                //              cellLocZ + k1, 1);
                                indy2 = get_col_index_Lx(
                                    cellLocX05 - cellLocX + i1,
                                    cellLocY05 - cellLocY + j1, k1, Y);
                                LmatX2[indx1][indy2] += value;

                            // xz
                            value = wx * wz1 * Ap * alpha *
                                    (-h.y() + alpha * h.x() * h.z());

                                indx1 = vind(cellLocX05 + i, cellLocY + j,
                                             cellLocZ + k, 0);
                                // indz2 = vind(cellLocX + i1, cellLocY + j1,
                                //              cellLocZ05 + k1, 2);
                                indz2 = get_col_index_Lx(i1,
                                                         j1,
                                                         cellLocZ05 - cellLocZ + k1, Z);
                                LmatX[indx1][indz2] += value;

                            // yx
                            value = wy * wx1 * Ap * alpha *
                                    (-h.z() + alpha * h.x() * h.y());

                                indy1 = vind(cellLocX + i, cellLocY05 + j,
                                             cellLocZ + k, 1);
                                // indx2 = vind(cellLocX05 + i1, cellLocY + j1,
                                //              cellLocZ + k1, 0);
                                indz2 = get_col_index_Ly(
                                    i1 + cellLocX05 - cellLocX, j1, k1, X);
                                LmatX[indy1][indx2] += value;

                            // yy
                            value = wy * wy1 * Ap *
                                    (1. + alpha * alpha * h.y() * h.y());

                                indy1 = vind(cellLocX + i, cellLocY05 + j,
                                             cellLocZ + k, 1);
                                // indy2 = vind(cellLocX + i1, cellLocY05 + j1,
                                //              cellLocZ + k1, 1);
                                indz2 = get_col_index_Ly(
                                    i1 , j1, k1, Y);
                                LmatX[indy1][indy2] += value;

                            // yz
                            value = wy * wz1 * Ap * alpha *
                                    (h.x() + alpha * h.y() * h.z());

                                indy1 = vind(cellLocX + i, cellLocY05 + j,
                                             cellLocZ + k, 1);
                                indz2 = vind(cellLocX + i1, cellLocY + j1,
                                             cellLocZ05 + k1, 2);
                                indz2 = get_col_index_Ly(
                                    i1, j1, k1 + cellLocZ05 - cellLocZ, Z);
                                LmatX[indy1][indz2] += value;

                            // zx
                            value = wz * wx1 * Ap * alpha *
                                    (h.y() + alpha * h.x() * h.z());

                                indz1 = vind(cellLocX + i, cellLocY + j,
                                             cellLocZ05 + k, 2);
                                // indx2 = vind(cellLocX05 + i1, cellLocY + j1,
                                //              cellLocZ + k1, 0);
                                indz2 = get_col_index_Lz(
                                    i1 + cellLocX05 - cellLocX, j1, k1, X);
                                LmatX[indz1][indx2] += value;

                            value = wz * wy1 * Ap * alpha *
                                    (-h.x() + alpha * h.y() * h.z());
                                indz1 = vind(cellLocX + i, cellLocY + j,
                                             cellLocZ05 + k, 2);
                                // indy2 = vind(cellLocX + i1, cellLocY05 + j1,
                                //              cellLocZ + k1, 1);
                                indz2 = get_col_index_Lz(
                                    i1, j1 + cellLocY05 - cellLocY, k1, Y);
                                LmatX[indz1][indy2] += value;

                            value = wz * wz1 * Ap *
                                    (1. + alpha * alpha * h.z() * h.z());
                            //if (fabs(value) > 1.e-16) {
                                indz1 = vind(cellLocX + i, cellLocY + j,
                                             cellLocZ05 + k, 2);
                                // indz2 = vind(cellLocX + i1, cellLocY + j1,
                                //              cellLocZ05 + k1, 2);
                                indz2 = get_col_index_Lz(
                                    i1, j1, k1 + cellLocZ05 - cellLocZ, Z);
                                LmatX[indz1][indz2] += value;
                          //  }
                        }   // G'z
                    }       // G'y
                }           // G'x
            }               // Gz
        }                   // Gy
    }                       // Gx
}

void Mesh::apply_periodic_boundaries(std::vector<IndexMap>& LmatX) {
    const auto size = int3(xSize, ySize, zSize);
    constexpr int OVERLAP_SIZE = 3;
    const int last_indx = size.x() - OVERLAP_SIZE;
    const int last_indy = size.y() - OVERLAP_SIZE;
    const int last_indz = size.z() - OVERLAP_SIZE;

    if (bounds.isPeriodic(X)) {
#pragma omp parallel for
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);

            if (ix < OVERLAP_SIZE) {
                for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                    auto ind2 = it->first;
                    auto value = it->second;
                    auto ix1 = pos_vind(ind2, 0);
                    auto iy1 = pos_vind(ind2, 1);
                    auto iz1 = pos_vind(ind2, 2);
                    auto id1 = pos_vind(ind2, 3);
                    auto indBound = vind(last_indx + ix, iy, iz, id);

                    if (ix1 < OVERLAP_SIZE) {
                        auto indBound2 = vind(last_indx + ix1, iy1, iz1, id1);
                        LmatX[indBound][indBound2] += value;
                    } else {
                        LmatX[indBound][ind2] += value;
                    }
                }
            }
        }
    }

    if (bounds.isPeriodic(Y)) {
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( iy < OVERLAP_SIZE ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(ix,last_indy + iy,iz,id);
            if(iy1 < OVERLAP_SIZE){
              auto indBound2 = vind(ix1,last_indy + iy1,iz1,id1);
              LmatX[indBound][indBound2] += value;
            } else {
              LmatX[indBound][ind2] += value;              
            }
          }
        }
    }
    }

    if (bounds.isPeriodic(Z)) {
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( iz < OVERLAP_SIZE ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(ix,iy,last_indz + iz,id);
            if(iz1 < OVERLAP_SIZE){
              auto indBound2 = vind(ix1,iy1,last_indz + iz1,id1);
              LmatX[indBound][indBound2] += value;
            } else {
              LmatX[indBound][ind2] += value;              
            }
          }
        }
    }
    }

    if (bounds.isPeriodic(X)) {
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( ix > last_indx - 1 ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(ix - last_indx,iy,iz,id);            
            
            if(ix1 > last_indx - 1){
              auto indBound2 = vind(ix1 - last_indx,iy1,iz1,id1);
              LmatX[indBound][indBound2] = value;
            } else{
              LmatX[indBound][ind2] = value;              
            }
          }
        }
    }
    }

    if (bounds.isPeriodic(Y)) {
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( iy > last_indy - 1 ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(ix, iy - last_indy,iz,id);            

            if(iy1 > last_indy - 1){
              auto indBound2 = vind(ix1, iy1 - last_indy,iz1,id1);
              LmatX[indBound][indBound2] = value;
            } else{
              LmatX[indBound][ind2] = value;              
            }
          }
        }      
    }
    }

    if (bounds.isPeriodic(Z)) {
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( iz > last_indz - 1 ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(ix, iy ,iz- last_indz,id);            

            if(iz1 > last_indz - 1){
              auto indBound2 = vind(ix1, iy1,iz1- last_indz,id1);
              LmatX[indBound][indBound2] = value;
            } else{
              LmatX[indBound][ind2] = value;              
            }
          }
        }      
    }
    }
}

void apply_periodic_border_with_add(Field3d &field, const Bounds &bounds) {

  auto size = field.size();
  auto nd = field.nd();

  constexpr int OVERLAP_SIZE = 3;
  auto i_max = size.x() - OVERLAP_SIZE;
  auto j_max = size.y() - OVERLAP_SIZE; 
  auto k_max = size.z() - OVERLAP_SIZE;

  if (bounds.isPeriodic(X)) {
      for (auto i = 0; i < OVERLAP_SIZE; ++i) {
          for (auto j = 0; j < size.y(); ++j) {
              for (auto k = 0; k < size.z(); ++k) {
                  for (auto dim = 0; dim < nd; dim++) {
                      field(i, j, k, dim) += field(i + i_max, j, k, dim);
                      field(i + i_max, j, k, dim) = field(i, j, k, dim);
                  }
              }
          }
      }
  }

  if (bounds.isPeriodic(Y)) {
      for (auto i = 0; i < size.x(); ++i) {
          for (auto j = 0; j < OVERLAP_SIZE; ++j) {
              for (auto k = 0; k < size.z(); ++k) {
                  for (auto dim = 0; dim < nd; dim++) {
                      field(i, j, k, dim) += field(i, j + j_max, k, dim);
                      field(i, j + j_max, k, dim) = field(i, j, k, dim);
                  }
              }
          }
      }
  }
  if (bounds.isPeriodic(Z)) {
      for (auto i = 0; i < size.x(); ++i) {
          for (auto j = 0; j < size.y(); ++j) {
              for (auto k = 0; k < OVERLAP_SIZE; ++k) {
                  for (auto dim = 0; dim < nd; dim++) {
                      field(i, j, k, dim) += field(i, j, k + k_max, dim);
                      field(i, j, k + k_max, dim) = field(i, j, k, dim);
                  }
              }
          }
      }
  }
}

void Mesh::apply_periodic_boundaries(Field3d& field) {
    apply_periodic_border_with_add(field, bounds);
}

void Mesh::apply_open_boundaries_z(Field3d& field) {
    constexpr int BOUNDARY_MARGIN = 2;
    auto size = field.size();
    auto nd = field.nd();

    if (bounds.lowerBounds.z == BoundType::OPEN) {
        for (auto ix = 0; ix < size.x(); ++ix) {
            for (auto iy = 0; iy < size.y(); ++iy) {
                field(ix, iy, 0, X) = 0.;
                field(ix, iy, 1, X) = 0.;

                field(ix, iy, 0, Y) = 0.;
                field(ix, iy, 1, Y) = 0.;

                field(ix, iy, 0, Z) = 0.;
            }
        }
    }
    if (bounds.upperBounds.z == BoundType::OPEN) {
        for (auto ix = 0; ix < size.x(); ++ix) {
            for (auto iy = 0; iy < size.y(); ++iy) {
                for (auto iz = size.z() - BOUNDARY_MARGIN; iz < size.z();
                     ++iz) {
                    for (auto dim = 0; dim < nd; dim++) {
                        field(ix, iy, iz, dim) = 0.;
                    }
                }
            }
        }
    }
}


void Mesh::apply_boundaries(Field3d& field) {
    apply_periodic_boundaries(field);
    apply_open_boundaries_z(field);
}

void Mesh::apply_boundaries(std::vector<IndexMap>& LmatX) {
    apply_periodic_boundaries(LmatX);
    apply_open_boundaries_z(LmatX);
}

void Mesh::apply_open_boundaries_z(std::vector<IndexMap>& LmatX) {
    constexpr int BOUNDARY_MARGIN = 2;

    const auto size = int3(xSize, ySize, zSize);
#pragma omp parallel for
    for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
        auto iz = pos_vind(i, Z);
        auto id = pos_vind(i, C);

        for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
            auto ind2 = it->first;
            auto iz1 = pos_vind(ind2, Z);
            auto id1 = pos_vind(ind2, C);
            if (bounds.lowerBounds.z == BoundType::OPEN) {
                if (iz < BOUNDARY_MARGIN || iz1 < BOUNDARY_MARGIN) {
                    if (!((iz == 1 && id == Z) || (iz1 == 1 && id1 == Z))) {
                        it->second = 0.;
                    }
                }
            }
            if (bounds.upperBounds.z == BoundType::OPEN) {
                if (iz >= size.z() - BOUNDARY_MARGIN ||
                    iz1 >= size.z() - BOUNDARY_MARGIN) {
                    it->second = 0.;
                }
            }
        }
    }
}

double3 Mesh::get_fieldE_in_cell(const Field3d& fieldE, int i,
                                                 int j, int k) const {
                    double3 E;
                    E.x() =
                        0.25 *
                        (fieldE(i, j, k, 0) + fieldE(i, j, k + 1, 0) +
                         fieldE(i, j + 1, k, 0) + fieldE(i, j + 1, k + 1, 0));

                    E.y() =
                        0.25 *
                        (fieldE(i, j, k, 1) + fieldE(i, j, k + 1, 1) +
                         fieldE(i + 1, j, k, 1) + fieldE(i + 1, j, k + 1, 1));

                    E.z() =
                        0.125 *
                        (fieldE(i, j, k, 2) + fieldE(i, j, k + 1, 2) +
                         fieldE(i, j + 1, k, 2) + fieldE(i, j + 1, k + 1, 2) +
                         fieldE(i + 1, j, k, 2) + fieldE(i + 1, j, k + 1, 2) +
                         fieldE(i + 1, j + 1, k, 2) +
                         fieldE(i + 1, j + 1, k + 1, 2));
                    return E;
}

double3 Mesh::get_fieldB_in_cell(const Field3d& fieldB, int i, int j,
                                 int k) const {
    double3 B;
    B.x() = 0.25 * (fieldB(i, j, k, 0) + fieldB(i, j, k + 1, 0) +
                    fieldB(i + 1, j, k, 0) + fieldB(i + 1, j, k + 1, 0));

    B.y() = 0.25 * (fieldB(i, j, k, 1) + fieldB(i, j, k + 1, 1) +
                    fieldB(i, j + 1, k, 1) + fieldB(i, j + 1, k + 1, 1));

    B.z() = fieldB(i, j, k, 2);
    return B;
}

double3 get_fieldE_in_pos(const Field3d& fieldE, const double3& coord, const Domain &domain) {
    double3 E;
    int indx, indy,indz,indx1, indy1,indz1;
    double sx0,sy0,sz0,sdx0,sdy0,sdz0;
    double sx1,sy1,sz1,sdx1,sdy1,sdz1;

    const double xx = coord.x() / domain.cell_size().x();
    const double yy = coord.y() / domain.cell_size().y();
    const double zz = coord.z() / domain.cell_size().z();

    indx = int(xx + 1.);
    indy = int(yy + 1.);
    indz = int(zz + 1.);


    indx1 = int(xx + 0.5);
    indy1 = int(yy + 0.5);
    indz1 = int(zz + 0.5); 
    
    sx1 = (xx - indx + 1.);
    sy1 = (yy - indy + 1.);
    sz1 = (zz - indz + 1.);
    sdx1 = (xx - indx1 + 0.5);
    sdy1 = (yy - indy1 + 0.5);
    sdz1 = (zz - indz1 + 0.5);


    sx0 = 1. - sx1;
    sy0 = 1. - sy1;
    sz0 = 1. - sz1;
    sdx0 = 1. - sdx1;
    sdy0 = 1. - sdy1;
    sdz0 = 1. - sdz1;

    E.x() = sdx0 * ( sy0 * ( sz0 * fieldE(indx1,indy,indz,0) + sz1 * fieldE(indx1,indy,indz+1,0) ) 
                + sy1 * ( sz0 * fieldE(indx1,indy+1,indz,0) + sz1 * fieldE(indx1,indy+1,indz+1,0) ) ) 
        + sdx1 * ( sy0 * ( sz0 * fieldE(indx1+1,indy,indz,0) + sz1 * fieldE(indx1+1,indy,indz+1,0) ) 
                + sy1 * ( sz0 * fieldE(indx1+1,indy+1,indz,0) + sz1 * fieldE(indx1+1,indy+1,indz+1,0) ) );

    E.y() = sx0 * ( sdy0 * ( sz0 * fieldE(indx,indy1,indz,1) + sz1 * fieldE(indx,indy1,indz+1,1) ) 
                + sdy1 * ( sz0 * fieldE(indx,indy1+1,indz,1) + sz1 * fieldE(indx,indy1+1,indz+1,1) ) ) 
        + sx1 * ( sdy0 * ( sz0 * fieldE(indx+1,indy1,indz,1) + sz1 * fieldE(indx+1,indy1,indz+1,1) ) 
                + sdy1 * ( sz0 * fieldE(indx+1,indy1+1,indz,1) + sz1 * fieldE(indx+1,indy1+1,indz+1,1) ) );

    E.z() = sx0 * ( sy0 * ( sdz0 * fieldE(indx,indy,indz1,2) + sdz1 * fieldE(indx,indy,indz1+1,2) ) 
                + sy1 * ( sdz0 * fieldE(indx,indy+1,indz1,2) + sdz1 * fieldE(indx,indy+1,indz1+1,2) ) ) 
        + sx1 * ( sy0 * ( sdz0 * fieldE(indx+1,indy,indz1,2) + sdz1 * fieldE(indx+1,indy,indz1+1,2) ) 
                + sy1 * ( sdz0 * fieldE(indx+1,indy+1,indz1,2) + sdz1 * fieldE(indx+1,indy+1,indz1+1,2) ) );

    return E;
}

double3 get_fieldB_in_pos(const Field3d& fieldB, const double3& coord, const Domain &domain) {

    double3 B;
    int indx, indy,indz,indx1, indy1,indz1;
    double sx0,sy0,sz0,sdx0,sdy0,sdz0;
    double sx1,sy1,sz1,sdx1,sdy1,sdz1;
 
    const double xx = coord.x() / domain.cell_size().x();
    const double yy = coord.y() / domain.cell_size().y();
    const double zz = coord.z() / domain.cell_size().z();
    
    indx = int(xx + 1.);
    indy = int(yy + 1.);
    indz = int(zz + 1.);


    indx1 = int(xx + 0.5);
    indy1 = int(yy + 0.5);
    indz1 = int(zz + 0.5); 
    
    sx1 = (xx - indx + 1.);
    sy1 = (yy - indy + 1.);
    sz1 = (zz - indz + 1.);
    sdx1 = (xx - indx1 + 0.5);
    sdy1 = (yy - indy1 + 0.5);
    sdz1 = (zz - indz1 + 0.5);


    sx0 = 1. - sx1;
    sy0 = 1. - sy1;
    sz0 = 1. - sz1;
    sdx0 = 1. - sdx1;
    sdy0 = 1. - sdy1;
    sdz0 = 1. - sdz1;

  B.x() = sx0 * ( sdy0 * ( sdz0 * fieldB(indx,indy1,indz1,0) + sdz1 * fieldB(indx,indy1,indz1+1,0) ) 
               + sdy1 * ( sdz0 * fieldB(indx,indy1+1,indz1,0) + sdz1 * fieldB(indx,indy1+1,indz1+1,0) ) ) 
       + sx1 * ( sdy0 * ( sdz0 * fieldB(indx+1,indy1,indz1,0) + sdz1 * fieldB(indx+1,indy1,indz1+1,0) ) 
               + sdy1 * ( sdz0 * fieldB(indx+1,indy1+1,indz1,0) + sdz1 * fieldB(indx+1,indy1+1,indz1+1,0) ) );

  B.y() = sdx0 * ( sy0 * ( sdz0 * fieldB(indx1,indy,indz1,1) + sdz1 * fieldB(indx1,indy,indz1+1,1) ) 
              + sy1 * ( sdz0 * fieldB(indx1,indy+1,indz1,1) + sdz1 * fieldB(indx1,indy+1,indz1+1,1) ) ) 
      + sdx1 * ( sy0 * ( sdz0 * fieldB(indx1+1,indy,indz1,1) + sdz1 * fieldB(indx1+1,indy,indz1+1,1) ) 
              + sy1 * ( sdz0 * fieldB(indx1+1,indy+1,indz1,1) + sdz1 * fieldB(indx1+1,indy+1,indz1+1,1) ) );

  B.z() = sdx0 * ( sdy0 * ( sz0 * fieldB(indx1,indy1,indz,2) + sz1 * fieldB(indx1,indy1,indz+1,2) ) 
              + sdy1 * ( sz0 * fieldB(indx1,indy1+1,indz,2) + sz1 * fieldB(indx1,indy1+1,indz+1,2) ) ) 
      + sdx1 * ( sdy0 * ( sz0 * fieldB(indx1+1,indy1,indz,2) + sz1 * fieldB(indx1+1,indy1,indz+1,2) ) 
              + sdy1 * ( sz0 * fieldB(indx1+1,indy1+1,indz,2) + sz1 * fieldB(indx1+1,indy1+1,indz+1,2) ) );

    return B;
}

// void get_fields_in_pos(const Field3d& fieldE,const Field3d& fieldB, const double3& r, double3& locE, double3 &locB) {

//     int indx, indy,indz,indx1, indy1,indz1;

//     double xx,yy,zz;

//     double sx0,sy0,sz0,sdx0,sdy0,sdz0;
//     double sx1,sy1,sz1,sdx1,sdy1,sdz1;

//     const double rdx = 1. / Dx;
//     const double rdy = 1. / Dy;
//     const double rdz = 1. / Dz;

//     xx = r.x() * rdx;
//     yy = r.y() * rdy;
//     zz = r.z() * rdz;
    
//     indx = int(xx + 1.);
//     indy = int(yy + 1.);
//     indz = int(zz + 1.);


//     indx1 = int(xx + 0.5);
//     indy1 = int(yy + 0.5);
//     indz1 = int(zz + 0.5); 
    
//     sx1 = (xx - indx + 1.);
//     sy1 = (yy - indy + 1.);
//     sz1 = (zz - indz + 1.);
//     sdx1 = (xx - indx1 + 0.5);
//     sdy1 = (yy - indy1 + 0.5);
//     sdz1 = (zz - indz1 + 0.5);


//     sx0 = 1. - sx1;
//     sy0 = 1. - sy1;
//     sz0 = 1. - sz1;
//     sdx0 = 1. - sdx1;
//     sdy0 = 1. - sdy1;
//     sdz0 = 1. - sdz1;
//     locE.x() = sdx0 * ( sy0 * ( sz0 * fieldE(indx1,indy,indz,0) + sz1 * fieldE(indx1,indy,indz+1,0) ) 
//                 + sy1 * ( sz0 * fieldE(indx1,indy+1,indz,0) + sz1 * fieldE(indx1,indy+1,indz+1,0) ) ) 
//         + sdx1 * ( sy0 * ( sz0 * fieldE(indx1+1,indy,indz,0) + sz1 * fieldE(indx1+1,indy,indz+1,0) ) 
//                 + sy1 * ( sz0 * fieldE(indx1+1,indy+1,indz,0) + sz1 * fieldE(indx1+1,indy+1,indz+1,0) ) );

//     locE.y() = sx0 * ( sdy0 * ( sz0 * fieldE(indx,indy1,indz,1) + sz1 * fieldE(indx,indy1,indz+1,1) ) 
//                 + sdy1 * ( sz0 * fieldE(indx,indy1+1,indz,1) + sz1 * fieldE(indx,indy1+1,indz+1,1) ) ) 
//         + sx1 * ( sdy0 * ( sz0 * fieldE(indx+1,indy1,indz,1) + sz1 * fieldE(indx+1,indy1,indz+1,1) ) 
//                 + sdy1 * ( sz0 * fieldE(indx+1,indy1+1,indz,1) + sz1 * fieldE(indx+1,indy1+1,indz+1,1) ) );

//     locE.z() = sx0 * ( sy0 * ( sdz0 * fieldE(indx,indy,indz1,2) + sdz1 * fieldE(indx,indy,indz1+1,2) ) 
//                 + sy1 * ( sdz0 * fieldE(indx,indy+1,indz1,2) + sdz1 * fieldE(indx,indy+1,indz1+1,2) ) ) 
//         + sx1 * ( sy0 * ( sdz0 * fieldE(indx+1,indy,indz1,2) + sdz1 * fieldE(indx+1,indy,indz1+1,2) ) 
//                 + sy1 * ( sdz0 * fieldE(indx+1,indy+1,indz1,2) + sdz1 * fieldE(indx+1,indy+1,indz1+1,2) ) );

//   locB.x() = sx0 * ( sdy0 * ( sdz0 * fieldB(indx,indy1,indz1,0) + sdz1 * fieldB(indx,indy1,indz1+1,0) ) 
//                + sdy1 * ( sdz0 * fieldB(indx,indy1+1,indz1,0) + sdz1 * fieldB(indx,indy1+1,indz1+1,0) ) ) 
//        + sx1 * ( sdy0 * ( sdz0 * fieldB(indx+1,indy1,indz1,0) + sdz1 * fieldB(indx+1,indy1,indz1+1,0) ) 
//                + sdy1 * ( sdz0 * fieldB(indx+1,indy1+1,indz1,0) + sdz1 * fieldB(indx+1,indy1+1,indz1+1,0) ) );

//   locB.y() = sdx0 * ( sy0 * ( sdz0 * fieldB(indx1,indy,indz1,1) + sdz1 * fieldB(indx1,indy,indz1+1,1) ) 
//               + sy1 * ( sdz0 * fieldB(indx1,indy+1,indz1,1) + sdz1 * fieldB(indx1,indy+1,indz1+1,1) ) ) 
//       + sdx1 * ( sy0 * ( sdz0 * fieldB(indx1+1,indy,indz1,1) + sdz1 * fieldB(indx1+1,indy,indz1+1,1) ) 
//               + sy1 * ( sdz0 * fieldB(indx1+1,indy+1,indz1,1) + sdz1 * fieldB(indx1+1,indy+1,indz1+1,1) ) );

//   locB.z() = sdx0 * ( sdy0 * ( sz0 * fieldB(indx1,indy1,indz,2) + sz1 * fieldB(indx1,indy1,indz+1,2) ) 
//               + sdy1 * ( sz0 * fieldB(indx1,indy1+1,indz,2) + sz1 * fieldB(indx1,indy1+1,indz+1,2) ) ) 
//       + sdx1 * ( sdy0 * ( sz0 * fieldB(indx1+1,indy1,indz,2) + sz1 * fieldB(indx1+1,indy1,indz+1,2) ) 
//               + sdy1 * ( sz0 * fieldB(indx1+1,indy1+1,indz,2) + sz1 * fieldB(indx1+1,indy1+1,indz+1,2) ) );

// }
