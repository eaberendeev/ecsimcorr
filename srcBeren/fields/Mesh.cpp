#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "solverSLE.h"
#include "util.h"

void Mesh::init(const Domain &domain, const ParametersMap &parameters){
    bounds.setBounds(domain.lower_bounds(), domain.upper_bounds());

    Lmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Lmat2.resize(domain.total_size() * 3, domain.total_size() * 3);
    Mmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Imat.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlE.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlB.resize(domain.total_size() * 3, domain.total_size() * 3);
    IMmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    LmatX.resize(domain.total_size() * 3);
    LmatX2.resize(domain.total_size() * 3);
    chargeDensityOld.resize(domain.size(), 1);
    chargeDensity.resize(domain.size(), 1);
    divE.resize(domain.total_size(), domain.total_size() * 3);
    LmatX2.resize(domain.total_size());
    // LmatX2.reserve();

    xCellSize = domain.cell_size().x();
    yCellSize = domain.cell_size().y();
    zCellSize = domain.cell_size().z();
    xSize = domain.size().x();
    ySize = domain.size().y();
    zSize = domain.size().z();

    stencil_Imat(Imat, domain);
    stencil_curlE(curlE, domain);
    stencil_curlB(curlB, domain);

    Operator curlB2(domain.total_size() * 3, domain.total_size() * 3);
    Operator curlE2(domain.total_size() * 3, domain.total_size() * 3);
    stencil_curlB_openZ(curlB2, domain);
   stencil_curlE_openZ(curlE2, domain);
   std::cout <<" norms curl " << (curlE-curlE2).norm() <<" " << (curlB - curlB2).norm() <<  std::endl;

    stencil_divE(divE, domain);

    double dt = parameters.get_double("Dt");
    Mmat = -0.25 * dt * dt * curlB * curlE;
    IMmat = Imat - Mmat;
    IMmat.makeCompressed();
}

void Mesh::zeroBoundL(Operator & mat) {
    const auto size = int3(xSize, ySize, zSize);

    mat.makeCompressed();
    // Получаем указатели на внутренние данные CSR
    double* values = mat.valuePtr();          // Массив значений
    const int* outerIndex = mat.outerIndexPtr(); // Массив индексов начала строк
    const int* innerIndex = mat.innerIndexPtr(); // Массив индексов столбцов

// #pragma omp parallel for schedule(dynamic, 32)
//     for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//         auto iz = pos_vind(i, 2);
//         auto id = pos_vind(i, 3);

//         for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
//             auto ind2 = it->first;
//             auto iz1 = pos_vind(ind2, 2);
//             auto id1 = pos_vind(ind2, 3);
//             if (iz < 2 || iz1 < 2 || iz > size.z() - 3 || iz1 > size.z() - 3) {
//                 if (!((iz == 1 && id == 2) || (iz1 == 1 && id1 == 2))) {
//                     it->second = 0.;
//                 }
//             }
//         }
//     }

#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
        auto iz = pos_vind(i, 2);
        auto id = pos_vind(i, 3);
        // Получаем диапазон ненулевых элементов для строки i
        int rowStart = outerIndex[i];
        int rowEnd = outerIndex[i + 1];

        for (int j = rowStart; j < rowEnd; j++) {

            int col = innerIndex[j];   // Столбец текущего элемента
            auto ind2 = col;
            auto iz1 = pos_vind(ind2, 2);
            auto id1 = pos_vind(ind2, 3);
            if (iz < 2 || iz1 < 2 || iz > size.z() - 3 || iz1 > size.z() - 3) {
                if (!((iz == 1 && id == 2) || (iz1 == 1 && id1 == 2))) {
                    values[j] = 0;
                }
            }
        }
    }
}

void Mesh::zeroBoundJ(Field3d& field)
{
      const auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32)
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
            if(iz < 2 || iz > size.z()-3) {
              if( !( (iz == 1 && id==2) )){
                field(ix,iy,iz,id) = 0.;
              }
            }
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
    auto sizes = field.sizes();
    for (auto i = 0; i < sizes.x(); ++i) {
        for (auto j = 0; j < sizes.y(); ++j) {
            for (auto k = 0; k < sizes.z(); ++k) {
                field(i, j, k, Dim::X) = bx;
                field(i, j, k, Dim::Y) = by;
                field(i, j, k, Dim::Z) = bz;
            }
        }
    }
}

void Mesh::prepare(){
#pragma omp parallel for
  for ( size_t i = 0; i < LmatX.size(); i++){
      for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
        it->second = 0.;
    }
  }
}

double Mesh::calc_energy_field(const Field3d& field) const{
  double potE = 0;
  const auto sizes = field.sizes();
  int i_max = sizes.x();
  int j_max = sizes.y();
  int k_max = sizes.z();
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
  const auto sizes = fieldE.sizes();
  int i_max = sizes.x();
  int j_max = sizes.y();
  int k_max = sizes.z();

  constexpr int OVERLAP_SIZE = 3;
  if (bounds.isPeriodic(X)) {
      i_max -= OVERLAP_SIZE;
  }
  if (bounds.isPeriodic(Y)) {
      j_max -= OVERLAP_SIZE;
  }
  if (false && bounds.isPeriodic(Z)) {
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
    const auto sizes = fieldE.sizes();
    int i_max = sizes.x();
    int j_max = sizes.y();
    int k_max = sizes.z();

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
                    Field3d& J, const double dt) {
    zeroBoundJ(J);

    Field3d rhs = E - dt*J + dt*curlB*B + Mmat*E;

    // solve Ax=b, fieldEn - output
    //double time1 = omp_get_wtime();
    //solve_linear_system<bicgstab>(IMmat, rhs.data(), En.data(), E.data());
    double time11 = omp_get_wtime();

    //solve_linear_system_mix(A, rhs, En.data(), E.data());

    solve_linear_system<BicgstabSolver<Field3d>>(IMmat, rhs, En, E);
    double time2 = omp_get_wtime();
    // solve_amgcl<Operator>(A, rhs, En.data(), E.data());
    //    double time2 = omp_get_wtime();
    std::cout<< "Correction fieldE solver error = "<< (Imat*En - Mmat*En - rhs).norm() << "\n";
    std::cout << "Correction fieldE Mysolver time = " << (time2 - time11) << "\n";
    //std::cout << "Correction fieldE Eigsolver time = " << (time11 - time1) << "\n";
    // std::cout<< "Correction fieldE solver time amgcl = "<< (time3-time2) <<
    // "\n";
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

    solve_linear_system<BicgstabSolver<Field>>(
        A, rhs, Enew.data(), E.data());

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
                    Field3d& J, double dt) {
    zeroBoundJ(J);
    zeroBoundL(Lmat2);

    double time1 = omp_get_wtime();

    Operator A = IMmat + Lmat2;
    double time11 = omp_get_wtime();
    Lmat2.makeCompressed();
 //   Operator A2 = parallel_sparse_addition2(IMmat, 1., Lmat, 1.);   // IMmat + Lmat;
  //  double time12 = omp_get_wtime();
   // std:: cout << "norm addition "<< (A-A2).norm() << "\n";
    //double time13 = omp_get_wtime();

    // std::cout << "Substraction matrix time = " << (omp_get_wtime() - time1)
    //           << "\n";

    //A.makeCompressed();

    // Field rhs = E.data() - dt * J.data() + dt * curlB * B.data() + L * E.data();
    // solve_linear_system<Operator, Field, BicgstabSolver<Operator, Field> >(A, rhs, Ep.data(),
    //                                                         E.data());
    //Field3d rhs = E - dt * J + dt * curlB * B + L * E;
    Field3d rhs = E - dt * J + dt * curlB * B + Mmat * E - Lmat2 * E;
    // auto matrix_op = [&](const Field3d& v) {
    //     return v - mv_product(Mmat, 1, Lmat, -1, v);
    // };
    // auto matrix_op = [](const Field3d& v, Field3d& res) {
    //     spmv(Mmat, 1, Lmat, -1, v, res);
    //     res = v - res;
    // };
    double time2 = omp_get_wtime();

    solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);

    double time21 = omp_get_wtime();
//     solve_amgcl<Operator>(A2, rhs.data(), Ep.data(), E.data());

    double time22 = omp_get_wtime();
    std::cout << "Prediction fieldE solver error = "
              << (IMmat  * Ep + Lmat2 * Ep - rhs).norm() << "\n";
     std::cout << "Prediction fieldE add matrices time = " << (time11 - time1) << "\n";
    //  std::cout << "Prediction fieldE parallel add matrices time = " << (time12 - time11)
    //            << "\n";
     std::cout << "Prediction fieldE Mysolver time = " << (time21 - time2)
               << "\n";
    //  std::cout << "Prediction fieldE amgcl time = " << (time22 - time21)
    //            << "\n";
     // std::cout << "Prediction fieldE solver time amgcl = " << (time3 - time2)
     //           << "\n";
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
       // B = double3(1,1,1);
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
                                    // std::cout << "LmatX[" << cellLocX05 + i << ", " << cellLocY + j << ", " << cellLocZ+k
                                    //           << "][" << cellLocX +i1 << ", " << cellLocY05 + j1 << ", " << cellLocZ+k1
                                    //           << "] += " << value << std::endl;
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

void Mesh::update_Lmat2(const double3& coord, const Domain &domain, double charge, double mass,
                       double mpw, const Field3d& fieldB, const double dt) {
    const int SMAX = 2; // SHAPE_SIZE;
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
       // B = double3(1,1,1);
        const double3 h = unit(B);

        const double alpha = 0.5*dt*charge*mag(B) / mass;
        const double Ap = 0.25*dt*dt*mpw*charge*charge / mass / (1+alpha*alpha);
        const int blockIndex = sind(cellLocX, cellLocY, cellLocZ);
        auto& currentBlock = LmatX2[blockIndex];

        const int xOffset = cellLocX05 - cellLocX + 1;
        const int yOffset = cellLocY05 - cellLocY + 1;
        const int zOffset = cellLocZ05 - cellLocZ + 1;
        for(i = 0; i < SMAX; ++i){
            for(j = 0; j < SMAX; ++j){
                for(k = 0; k < SMAX; ++k){
                    wx = sx05[i] * sy[j] * sz[k];
                    wy = sx[i] * sy05[j] * sz[k];
                    wz = sx[i] * sy[j] * sz05[k];
                    indx1 = indX(xOffset + i, j, k);
                    indy1 = indY(i, yOffset + j, k);
                    indz1 = indZ(i, j, zOffset + k);
                    for (i1 = 0; i1 < SMAX; ++i1) {
                        for (j1 = 0; j1 < SMAX; ++j1) {
                            for (k1 = 0; k1 < SMAX; ++k1) {
                                wx1 = sx05[i1] * sy[j1] * sz[k1];
                                wy1 = sx[i1] * sy05[j1] * sz[k1];
                                wz1 = sx[i1] * sy[j1] * sz05[k1];
                                indx2 = indX(xOffset + i1, j1, k1);
                                indy2 = indY(i1, yOffset + j1, k1);
                                indz2 = indZ(i1, j1, zOffset + k1);

                                // xx
                                value = wx * wx1 * Ap *
                                        (1. + alpha * alpha * h.x() * h.x());

                                    currentBlock(indx1, indx2, 0) += value;
                                    // xy
                                    value = wx * wy1 * Ap *
                                            (alpha * h.z() +
                                             alpha * alpha * h.x() * h.y());
                                    currentBlock(indx1, indy2, 1) += value;
                                    // xz
                                    value = wx * wz1 * Ap * alpha *
                                            (-h.y() + alpha * h.x() * h.z());
                                    currentBlock(indx1, indz2, 2) += value;
                                    // yx
                                    value = wy * wx1 * Ap * alpha *
                                            (-h.z() + alpha * h.x() * h.y());

                                    currentBlock(indy1, indx2, 3) += value;
                                    // yy
                                    value =
                                        wy * wy1 * Ap *
                                        (1. + alpha * alpha * h.y() * h.y());
                                    currentBlock(indy1, indy2, 4) += value;
                                    // yz
                                    value = wy * wz1 * Ap * alpha *
                                            (h.x() + alpha * h.y() * h.z());
                                    currentBlock(indy1, indz2, 5) += value;
                                    // zx
                                    value = wz * wx1 * Ap * alpha *
                                            (h.y() + alpha * h.x() * h.z());
                                    currentBlock(indz1, indx2, 6) += value;
                                    // zy
                                    value = wz * wy1 * Ap * alpha *
                                            (-h.x() + alpha * h.y() * h.z());
                                    currentBlock(indz1, indy2, 7) += value;
                                    value =
                                        wz * wz1 * Ap *
                                        (1. + alpha * alpha * h.z() * h.z());
                                    currentBlock(indz1, indz2, 8) += value;
                            }   // G'z
                        }   // G'y
                    }   // G'x
                } // Gz
            } // Gy
        } // Gx
}

void Mesh::update_LmatNGP(const double3& coord, const Domain& domain,
                       double charge, double mass, double mpw,
                       const Field3d& fieldB, const double dt) {
    double3 B;
    double value;
    int cellLocX, cellLocY, cellLocZ, cellLocX05, cellLocY05, cellLocZ05;
    double coordLocX, coordLocY, coordLocZ;
    double coordLocX05, coordLocY05, coordLocZ05;
    int indx1, indy1, indz1;
    int indx2, indy2, indz2;

    B = 0.;

    coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    coordLocX05 = coordLocX - 0.5;
    coordLocY05 = coordLocY - 0.5;
    coordLocZ05 = coordLocZ - 0.5;

    cellLocX = ngp(coordLocX);
    cellLocY = ngp(coordLocY);
    cellLocZ = ngp(coordLocZ);
    cellLocX05 = ngp(coordLocX05);
    cellLocY05 = ngp(coordLocY05);
    cellLocZ05 = ngp(coordLocZ05);

    B.x() = fieldB(cellLocX, cellLocY05, cellLocZ05, 0);
    B.y() = fieldB(cellLocX05, cellLocY, cellLocZ05, 1);
    B.z() = fieldB(cellLocX05, cellLocY05, cellLocZ, 2);

    // B = double3(1,1,1);
    const double3 h = unit(B);

    const double alpha = 0.5 * dt * charge * mag(B) / mass;
    const double Ap =
        0.25 * dt * dt * mpw * charge * charge / mass / (1 + alpha * alpha);

    value = Ap * (1. + alpha * alpha * h.x() * h.x());
    if (fabs(value) > 1.e-16) {
        indx1 = vind(cellLocX05, cellLocY, cellLocZ, 0);
        indx2 = vind(cellLocX05, cellLocY, cellLocZ, 0);
        LmatX[indx1][indx2] += value;
    }
    // xy
    value = Ap * (alpha * h.z() + alpha * alpha * h.x() * h.y());
    if (fabs(value) > 1.e-16) {
        indx1 = vind(cellLocX05, cellLocY, cellLocZ, 0);
        indy2 = vind(cellLocX, cellLocY05, cellLocZ, 1);
        LmatX[indx1][indy2] += value;
    }
    // xz
    value = Ap * alpha * (-h.y() + alpha * h.x() * h.z());
    if (fabs(value) > 1.e-16) {
        indx1 = vind(cellLocX05, cellLocY, cellLocZ, 0);
        indz2 = vind(cellLocX, cellLocY, cellLocZ05, 2);
        LmatX[indx1][indz2] += value;
    }
    // yx
    value = Ap * alpha * (-h.z() + alpha * h.x() * h.y());
    if (fabs(value) > 1.e-16) {
        indy1 = vind(cellLocX, cellLocY05, cellLocZ, 1);
        indx2 = vind(cellLocX05, cellLocY, cellLocZ, 0);
        LmatX[indy1][indx2] += value;
    }
    // yy
    value = Ap * (1. + alpha * alpha * h.y() * h.y());
    if (fabs(value) > 1.e-16) {
        indy1 = vind(cellLocX, cellLocY05, cellLocZ, 1);
        indy2 = vind(cellLocX, cellLocY05, cellLocZ, 1);
        LmatX[indy1][indy2] += value;
    }
    // yz
    value = Ap * alpha * (h.x() + alpha * h.y() * h.z());
    if (fabs(value) > 1.e-16) {
        indy1 = vind(cellLocX, cellLocY05, cellLocZ, 1);
        indz2 = vind(cellLocX, cellLocY, cellLocZ05, 2);
        LmatX[indy1][indz2] += value;
    }
    // zx
    value = Ap * alpha * (h.y() + alpha * h.x() * h.z());
    if (fabs(value) > 1.e-16) {
        indz1 = vind(cellLocX, cellLocY, cellLocZ05, 2);
        indx2 = vind(cellLocX05, cellLocY, cellLocZ, 0);
        LmatX[indz1][indx2] += value;
    }
    value = Ap * alpha * (-h.x() + alpha * h.y() * h.z());
    if (fabs(value) > 1.e-16) {
        indz1 = vind(cellLocX, cellLocY, cellLocZ05, 2);
        indy2 = vind(cellLocX, cellLocY05, cellLocZ, 1);
        LmatX[indz1][indy2] += value;
    }
    value = Ap * (1. + alpha * alpha * h.z() * h.z());
    if (fabs(value) > 1.e-16) {
        indz1 = vind(cellLocX, cellLocY, cellLocZ05, 2);
        indz2 = vind(cellLocX, cellLocY, cellLocZ05, 2);
        LmatX[indz1][indz2] += value;
    }
}

void Mesh::apply_periodic_boundaries(std::vector<IndexMap>& LmatX) {
    const auto size = int3(xSize, ySize, zSize);
    constexpr int OVERLAP_SIZE = 3;
    const int last_indx = size.x() - OVERLAP_SIZE;
    const int last_indy = size.y() - OVERLAP_SIZE;
    const int last_indz = size.z() - OVERLAP_SIZE;

    if (bounds.isPeriodic(X)) {
#pragma omp parallel for schedule(dynamic, 32)
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
#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);

            if (iy < OVERLAP_SIZE) {
                for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                    auto ind2 = it->first;
                    auto value = it->second;
                    auto ix1 = pos_vind(ind2, 0);
                    auto iy1 = pos_vind(ind2, 1);
                    auto iz1 = pos_vind(ind2, 2);
                    auto id1 = pos_vind(ind2, 3);
                    auto indBound = vind(ix, last_indy + iy, iz, id);
                    if (iy1 < OVERLAP_SIZE) {
                        auto indBound2 = vind(ix1, last_indy + iy1, iz1, id1);
                        LmatX[indBound][indBound2] += value;
                    } else {
                        LmatX[indBound][ind2] += value;
                    }
                }
            }
        }
    }

    if (bounds.isPeriodic(Z)) {
#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);

            if (iz < OVERLAP_SIZE) {
                for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                    auto ind2 = it->first;
                    auto value = it->second;
                    auto ix1 = pos_vind(ind2, 0);
                    auto iy1 = pos_vind(ind2, 1);
                    auto iz1 = pos_vind(ind2, 2);
                    auto id1 = pos_vind(ind2, 3);
                    auto indBound = vind(ix, iy, last_indz + iz, id);
                    if (iz1 < OVERLAP_SIZE) {
                        auto indBound2 = vind(ix1, iy1, last_indz + iz1, id1);
                        LmatX[indBound][indBound2] += value;
                    } else {
                        LmatX[indBound][ind2] += value;
                    }
                }
            }
        }
    }

    if (bounds.isPeriodic(X)) {
#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);

            if (ix > last_indx - 1) {
                for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                    auto ind2 = it->first;
                    auto value = it->second;
                    auto ix1 = pos_vind(ind2, 0);
                    auto iy1 = pos_vind(ind2, 1);
                    auto iz1 = pos_vind(ind2, 2);
                    auto id1 = pos_vind(ind2, 3);
                    auto indBound = vind(ix - last_indx, iy, iz, id);

                    if (ix1 > last_indx - 1) {
                        auto indBound2 = vind(ix1 - last_indx, iy1, iz1, id1);
                        LmatX[indBound][indBound2] = value;
                    } else {
                        LmatX[indBound][ind2] = value;
                    }
                }
            }
        }
    }

    if (bounds.isPeriodic(Y)) {
#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);

            if (iy > last_indy - 1) {
                for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                    auto ind2 = it->first;
                    auto value = it->second;
                    auto ix1 = pos_vind(ind2, 0);
                    auto iy1 = pos_vind(ind2, 1);
                    auto iz1 = pos_vind(ind2, 2);
                    auto id1 = pos_vind(ind2, 3);
                    auto indBound = vind(ix, iy - last_indy, iz, id);

                    if (iy1 > last_indy - 1) {
                        auto indBound2 = vind(ix1, iy1 - last_indy, iz1, id1);
                        LmatX[indBound][indBound2] = value;
                    } else {
                        LmatX[indBound][ind2] = value;
                    }
                }
            }
        }
    }

    if (bounds.isPeriodic(Z)) {
#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);

            if (iz > last_indz - 1) {
                for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                    auto ind2 = it->first;
                    auto value = it->second;
                    auto ix1 = pos_vind(ind2, 0);
                    auto iy1 = pos_vind(ind2, 1);
                    auto iz1 = pos_vind(ind2, 2);
                    auto id1 = pos_vind(ind2, 3);
                    auto indBound = vind(ix, iy, iz - last_indz, id);

                    if (iz1 > last_indz - 1) {
                        auto indBound2 = vind(ix1, iy1, iz1 - last_indz, id1);
                        LmatX[indBound][indBound2] = value;
                    } else {
                        LmatX[indBound][ind2] = value;
                    }
                }
            }
        }
    }
}

void apply_periodic_border_with_add(Field3d &field, const Bounds &bounds) {

  auto sizes = field.sizes();
  auto nd = field.nd();

  constexpr int OVERLAP_SIZE = 3;
  auto i_max = sizes.x() - OVERLAP_SIZE;
  auto j_max = sizes.y() - OVERLAP_SIZE; 
  auto k_max = sizes.z() - OVERLAP_SIZE;

  if (bounds.isPeriodic(X)) {
      for (auto i = 0; i < OVERLAP_SIZE; ++i) {
          for (auto j = 0; j < sizes.y(); ++j) {
              for (auto k = 0; k < sizes.z(); ++k) {
                  for (auto dim = 0; dim < nd; dim++) {
                      field(i, j, k, dim) += field(i + i_max, j, k, dim);
                      field(i + i_max, j, k, dim) = field(i, j, k, dim);
                  }
              }
          }
      }
  }

  if (bounds.isPeriodic(Y)) {
      for (auto i = 0; i < sizes.x(); ++i) {
          for (auto j = 0; j < OVERLAP_SIZE; ++j) {
              for (auto k = 0; k < sizes.z(); ++k) {
                  for (auto dim = 0; dim < nd; dim++) {
                      field(i, j, k, dim) += field(i, j + j_max, k, dim);
                      field(i, j + j_max, k, dim) = field(i, j, k, dim);
                  }
              }
          }
      }
  }
  if (bounds.isPeriodic(Z)) {
      for (auto i = 0; i < sizes.x(); ++i) {
          for (auto j = 0; j < sizes.y(); ++j) {
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

void Mesh::apply_open_boundaries(Field3d& field, const Domain& domain) {
    //    constexpr int BOUNDARY_MARGIN = 2;
    //    auto sizes = field.sizes();
    //auto nd = field.nd();
    auto size = field.sizes();
    auto setValuesZero =
        [](double& value, const Domain& domain, int vindg) {
            bool setZero = !domain.in_region_electric(vindg);
            if (setZero) {
                value = 0.;
            }
        };
#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < 3*(size.x() * size.y() * size.z()); i++) {
        setValuesZero(field[i], domain, i);
    }

    // if (bounds.lowerBounds.z == BoundType::OPEN) {
    //     for (auto ix = 0; ix < sizes.x(); ++ix) {
    //         for (auto iy = 0; iy < sizes.y(); ++iy) {
    //             field(ix, iy, 0, X) = 0.;
    //             field(ix, iy, 1, X) = 0.;

    //             field(ix, iy, 0, Y) = 0.;
    //             field(ix, iy, 1, Y) = 0.;

    //             field(ix, iy, 0, Z) = 0.;
    //         }
    //     }
    // }
    // if (bounds.upperBounds.z == BoundType::OPEN) {
    //     for (auto ix = 0; ix < sizes.x(); ++ix) {
    //         for (auto iy = 0; iy < sizes.y(); ++iy) {
    //             for (auto iz = sizes.z() - BOUNDARY_MARGIN; iz < sizes.z();
    //                  ++iz) {
    //                 for (auto dim = 0; dim < nd; dim++) {
    //                     field(ix, iy, iz, dim) = 0.;
    //                 }
    //             }
    //         }
    //     }
    // }
}

void Mesh::apply_boundaries(Field3d& field, const Domain& domain) {
    apply_periodic_boundaries(field);
    apply_open_boundaries(field, domain);
}

void Mesh::apply_density_periodic_boundaries(Field3d& field) {
    apply_periodic_border_with_add(field, bounds);
}

void Mesh::apply_density_open_boundaries(Field3d& field, const Domain& domain) {
//    constexpr int BOUNDARY_MARGIN = 2;
    auto size = field.sizes();
    auto setValuesZero =
        [](double& value, const Domain& domain, int vindg) {
            bool setZero = !domain.in_region_density(vindg);
            if (setZero) {
                value = 0.;
            }
        };
#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < (size.x() * size.y() * size.z()); i++) {
        setValuesZero(field[i], domain, i);
    }
    // if (bounds.lowerBounds.z == BoundType::OPEN) {
    //     for (auto ix = 0; ix < sizes.x(); ++ix) {
    //         for (auto iy = 0; iy < sizes.y(); ++iy) {
    //             field(ix, iy, 0, 0) = 0.;
    //             // field(ix, iy, 1, X) = 0.;
    //         }
    //     }
    // }
    // if (bounds.upperBounds.z == BoundType::OPEN) {
    //     for (auto ix = 0; ix < sizes.x(); ++ix) {
    //         for (auto iy = 0; iy < sizes.y(); ++iy) {
    //                 field(ix, iy, sizes.z() - 1, 0) = 0.;
    //         }
    //     }
    // }
}

void Mesh::apply_density_boundaries(Field3d& field, const Domain& domain) {
    apply_density_periodic_boundaries(field);
    apply_density_open_boundaries(field, domain);
}

void Mesh::apply_boundaries(std::vector<IndexMap>& LmatX,
                            const Domain& domain) {
    apply_periodic_boundaries(LmatX);
    apply_open_boundaries(LmatX, domain);
}

void Mesh::apply_open_boundaries(std::vector<IndexMap>& LmatX,
                                   const Domain& domain) {
    //constexpr int BOUNDARY_MARGIN = 2;

    const auto size = int3(xSize, ySize, zSize);
        auto setValuesZero =
        [](double& value, const Domain& domain, int vindg,
           int vindg1) {
            bool setZero = !domain.in_region_electric(vindg)
                                || !domain.in_region_electric(vindg1);
            if (setZero) {
                value = 0.;
            }
        };

#pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
            if (!domain.in_region_electric(i))
                continue;

            //auto iz = pos_vind(i, Z);
            //auto id = pos_vind(i, C);

            for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                auto ind2 = it->first;
                setValuesZero(it->second, domain, i, ind2);

                //auto iz1 = pos_vind(ind2, Z);
                //auto id1 = pos_vind(ind2, C);
                // if (bounds.lowerBounds.z == BoundType::OPEN) {
                //     if (iz < BOUNDARY_MARGIN || iz1 < BOUNDARY_MARGIN) {
                //         if (!((iz == 1 && id == Z) || (iz1 == 1 && id1 == Z))) {
                //             it->second = 0.;
                //         }
                //     }
                // }
                // if (bounds.upperBounds.z == BoundType::OPEN) {
                //     if (iz >= size.z() - BOUNDARY_MARGIN ||
                //         iz1 >= size.z() - BOUNDARY_MARGIN) {
                //         it->second = 0.;
                //     }
                // }
        }
    }
}

void Mesh::apply_open_boundaries_z(Field3d& field) {
    constexpr int BOUNDARY_MARGIN = 2;
    auto size = field.sizes();
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


void Mesh::apply_open_boundaries_z(std::vector<IndexMap>& mat) {
    constexpr int BOUNDARY_MARGIN = 2;

    const auto size = int3(xSize, ySize, zSize);
#pragma omp parallel for
    for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
        auto iz = pos_vind(i, Z);
        auto id = pos_vind(i, C);

        for (auto it = mat[i].begin(); it != mat[i].end(); ++it) {
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

// void Mesh::apply_open_boundaries_z(Operator& LmatX) {
//     constexpr int BOUNDARY_MARGIN = 2;

//     const auto size = int3(xSize, ySize, zSize);
// #pragma omp parallel for schedule(dynamic, 128)
//     for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//         auto iz = pos_vind(i, Z);
//         auto id = pos_vind(i, C);

//         for (Operator::InnerIterator it(LmatX, i); it; ++it) {
//             auto ind2 = it.col();
//             auto iz1 = pos_vind(ind2, Z);
//             auto id1 = pos_vind(ind2, C);
//             if (bounds.lowerBounds.z == BoundType::OPEN) {
//                 if (iz < BOUNDARY_MARGIN || iz1 < BOUNDARY_MARGIN) {
//                     if (!((iz == 1 && id == Z) || (iz1 == 1 && id1 == Z))) {
//                         LmatX.coeffRef(i, ind2) = 0.;
//                     }
//                 }
//             }
//             if (bounds.upperBounds.z == BoundType::OPEN) {
//                 if (iz >= size.z() - BOUNDARY_MARGIN ||
//                     iz1 >= size.z() - BOUNDARY_MARGIN) {
//                     LmatX.coeffRef(i, ind2) = 0.;
//                 }
//             }
//         }
//     }
// }

// void Mesh::apply_open_boundaries_z(Operator& LmatX) {
//     constexpr int BOUNDARY_MARGIN = 2;

//     const auto size = int3(xSize, ySize, zSize);
//     LmatX.makeCompressed();
//     // Получаем указатели на внутренние данные CSR
//     double* values = LmatX.valuePtr();          // Массив значений
//     const int* outerIndex = LmatX.outerIndexPtr(); // Массив индексов начала строк
//     const int* innerIndex = LmatX.innerIndexPtr(); // Массив индексов столбцов

// #pragma omp parallel for schedule(dynamic, 128)
//     for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//         auto iz = pos_vind(i, Z);
//         auto id = pos_vind(i, C);

//         // Получаем диапазон ненулевых элементов для строки i
//         int rowStart = outerIndex[i];
//         int rowEnd = outerIndex[i + 1];

//         for (int j = rowStart; j < rowEnd; j++) {
//             int col = innerIndex[j]; // Столбец текущего элемента
//             auto iz1 = pos_vind(col, Z);
//             auto id1 = pos_vind(col, C);

//             if (bounds.lowerBounds.z == BoundType::OPEN) {
//                 if (iz < BOUNDARY_MARGIN || iz1 < BOUNDARY_MARGIN) {
//                     if (!((iz == 1 && id == Z) || (iz1 == 1 && id1 == Z))) {
//                         values[j] = 0.; // Прямое изменение значения
//                     }
//                 }
//             }
//             if (bounds.upperBounds.z == BoundType::OPEN) {
//                 if (iz >= size.z() - BOUNDARY_MARGIN ||
//                     iz1 >= size.z() - BOUNDARY_MARGIN) {
//                     values[j] = 0.; // Прямое изменение значения
//                 }
//             }
//         }
//     }
// }

void Mesh::apply_open_boundaries(Operator& LmatX, Domain& domain) {

    const auto size = int3(xSize, ySize, zSize);
    LmatX.makeCompressed();
    // Получаем указатели на внутренние данные CSR
    double* values = LmatX.valuePtr();          // Массив значений
    const int* outerIndex = LmatX.outerIndexPtr(); // Массив индексов начала строк
    const int* innerIndex = LmatX.innerIndexPtr(); // Массив индексов столбцов

    auto setValuesZero =
        [](double* values, const Domain& domain, int vindg,
           int vindg1, int value_index) {
            bool setZero = !domain.in_region_electric(vindg)
                                || !domain.in_region_electric(vindg1);
            if (setZero) {
                values[value_index] = 0.;
            }
        };
#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
        if (!domain.in_region_electric(i)) continue;

        // Получаем диапазон ненулевых элементов для строки i
        int rowStart = outerIndex[i];
        int rowEnd = outerIndex[i + 1];

        for (int j = rowStart; j < rowEnd; j++) {
            int col = innerIndex[j]; // Столбец текущего элемента
            setValuesZero(values, domain, i, col, j);
        }
    }
}

void Mesh::apply_boundaries(Operator& LmatX, Domain& domain) {
    std::cout << "Applying periodic boundaries for CSR Lmat\n" << std::endl;
    apply_periodic_boundaries(LmatX);
    std::cout << "Applying open periodic boundaries for CSR Lmat\n" << std::endl;
    apply_open_boundaries(LmatX, domain);
    std::cout << "Boundary for CSR Lmat done!\n" << std::endl;
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

double3 interpolateE(const Field3d& fieldE, const double3& normalized_coord,
                     ShapeType type) {
    switch (type) {
        case ShapeType::NGP:
            return interpolateE_ngp(fieldE, normalized_coord);
        case ShapeType::Linear:
            return interpolateE_linear(fieldE, normalized_coord);
        case ShapeType::Quadratic:
            std::cout << "interpolateE for quadratic is not supported\n"
                      << std::endl;
    }
    return double3(0, 0, 0);
}

double3 interpolateE_ngp(const Field3d& fieldE, const double3& normalized_coord) {
    const auto x05 = normalized_coord.x() - 0.5;
    const auto y05 = normalized_coord.y() - 0.5;
    const auto z05 = normalized_coord.z() - 0.5;

    const auto ix = ngp(normalized_coord.x() + GHOST_CELLS);
    const auto iy = ngp(normalized_coord.y() + GHOST_CELLS);
    const auto iz = ngp(normalized_coord.z() + GHOST_CELLS);
    const auto ix05 = ngp(x05 + GHOST_CELLS);
    const auto iy05 = ngp(y05 + GHOST_CELLS);
    const auto iz05 = ngp(z05 + GHOST_CELLS);

    double3 E;
    E.x() = fieldE(ix05, iy, iz, X);
    E.y() = fieldE(ix, iy05, iz, Y);
    E.z() = fieldE(ix, iy, iz05, Z);
    return E;
}

double3 interpolateB_ngp(const Field3d& fieldB, const double3& normalized_coord) {
    const auto x05 = normalized_coord.x() - 0.5;
    const auto y05 = normalized_coord.y() - 0.5;
    const auto z05 = normalized_coord.z() - 0.5;

    const auto ix = ngp(normalized_coord.x() + GHOST_CELLS);
    const auto iy = ngp(normalized_coord.y() + GHOST_CELLS);
    const auto iz = ngp(normalized_coord.z() + GHOST_CELLS);
    const auto ix05 = ngp(x05 + GHOST_CELLS);
    const auto iy05 = ngp(y05 + GHOST_CELLS);
    const auto iz05 = ngp(z05 + GHOST_CELLS);

    double3 B;
    B.x() = fieldB(ix, iy05, iz05, X);
    B.y() = fieldB(ix05, iy, iz05, Y);
    B.z() = fieldB(ix05, iy05, iz, Z);
    return B;
}

double3 interpolateE_linear(const Field3d& fieldE, const double3& normalized_coord) {
    double3 E;

    const double xx = normalized_coord.x() + GHOST_CELLS;
    const double yy = normalized_coord.y() + GHOST_CELLS;
    const double zz = normalized_coord.z() + GHOST_CELLS;

    const int indx = int(xx);
    const int indy = int(yy);
    const int indz = int(zz);

    const int indx1 = int(xx - 0.5);
    const int indy1 = int(yy - 0.5);
    const int indz1 = int(zz - 0.5);

    const double sx1 = (xx - indx);
    const double sy1 = (yy - indy);
    const double sz1 = (zz - indz);
    const double sdx1 = (xx - indx1 - 0.5);
    const double sdy1 = (yy - indy1 - 0.5);
    const double sdz1 = (zz - indz1 - 0.5);

    const double sx0 = 1. - sx1;
    const double sy0 = 1. - sy1;
    const double sz0 = 1. - sz1;
    const double sdx0 = 1. - sdx1;
    const double sdy0 = 1. - sdy1;
    const double sdz0 = 1. - sdz1;

   E.x() = sdx0 * (sy0 * (sz0 * fieldE(indx1, indy, indz, 0) +
                          sz1 * fieldE(indx1, indy, indz + 1, 0)) +
                   sy1 * (sz0 * fieldE(indx1, indy + 1, indz, 0) +
                          sz1 * fieldE(indx1, indy + 1, indz + 1, 0))) +
           sdx1 * (sy0 * (sz0 * fieldE(indx1 + 1, indy, indz, 0) +
                          sz1 * fieldE(indx1 + 1, indy, indz + 1, 0)) +
                   sy1 * (sz0 * fieldE(indx1 + 1, indy + 1, indz, 0) +
                          sz1 * fieldE(indx1 + 1, indy + 1, indz + 1, 0)));

   E.y() = sx0 * (sdy0 * (sz0 * fieldE(indx, indy1, indz, 1) +
                          sz1 * fieldE(indx, indy1, indz + 1, 1)) +
                  sdy1 * (sz0 * fieldE(indx, indy1 + 1, indz, 1) +
                          sz1 * fieldE(indx, indy1 + 1, indz + 1, 1))) +
           sx1 * (sdy0 * (sz0 * fieldE(indx + 1, indy1, indz, 1) +
                          sz1 * fieldE(indx + 1, indy1, indz + 1, 1)) +
                  sdy1 * (sz0 * fieldE(indx + 1, indy1 + 1, indz, 1) +
                          sz1 * fieldE(indx + 1, indy1 + 1, indz + 1, 1)));

   E.z() = sx0 * (sy0 * (sdz0 * fieldE(indx, indy, indz1, 2) +
                         sdz1 * fieldE(indx, indy, indz1 + 1, 2)) +
                  sy1 * (sdz0 * fieldE(indx, indy + 1, indz1, 2) +
                         sdz1 * fieldE(indx, indy + 1, indz1 + 1, 2))) +
           sx1 * (sy0 * (sdz0 * fieldE(indx + 1, indy, indz1, 2) +
                         sdz1 * fieldE(indx + 1, indy, indz1 + 1, 2)) +
                  sy1 * (sdz0 * fieldE(indx + 1, indy + 1, indz1, 2) +
                         sdz1 * fieldE(indx + 1, indy + 1, indz1 + 1, 2)));

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
