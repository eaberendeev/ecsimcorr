#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "solverSLE.h"
#include "util.h"
#include "operators.h"
#include "interpolation.h"

void Mesh::init(const Domain &domain, double dt){
    bounds.setBounds(domain.lower_bounds(), domain.upper_bounds());

    Lmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Lmat2.resize(domain.total_size() * 3, domain.total_size() * 3);
    Mmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Imat.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlE.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlB.resize(domain.total_size() * 3, domain.total_size() * 3);
    IMmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    LmatX.resize(domain.total_size() * 3);
    chargeDensityOld.resize(domain.size(), 1);
    chargeDensity.resize(domain.size(), 1);
    divE.resize(domain.total_size(), domain.total_size() * 3);
    // TODO: move sind, converter func to BlockMatrix, resize with 3dim
    LmatX2.resize(domain.total_size());
    LmatX_NGP.resize(domain.total_size());
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

    Mmat = -0.25 * dt * dt * curlB * curlE;
    IMmat = Imat - Mmat;
    IMmat.makeCompressed();
}

void Mesh::zeroBoundL(Operator & mat) {
    const auto size = Vector3I(xSize, ySize, zSize);

    mat.makeCompressed();
    // Получаем указатели на внутренние данные CSR
    double* values = mat.valuePtr();          // Массив значений
    const int* outerIndex = mat.outerIndexPtr(); // Массив индексов начала строк
    const int* innerIndex = mat.innerIndexPtr(); // Массив индексов столбцов

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

void Mesh::zeroBoundJ(Field3d& field) {
    const auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
        auto ix = pos_vind(i, 0);
        auto iy = pos_vind(i, 1);
        auto iz = pos_vind(i, 2);
        auto id = pos_vind(i, 3);
        if (iz < 2 || iz > size.z() - 3) {
            if (!((iz == 1 && id == 2))) {
                field(ix, iy, iz, id) = 0.;
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

  for (auto i = 0; i < i_max; ++i) {
      for (auto j = 0; j < j_max; ++j) {
          for (auto k = 0; k < k_max; ++k) {
              Vector3R v = Vector3R(field(i, j, k, 0), field(i, j, k, 1),
                                  field(i, j, k, 2));
              potE += v.dot(v);
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
  if (bounds.isPeriodic(Z)) {
      k_max -= OVERLAP_SIZE;
  }
  for(auto i = 0; i < i_max; ++i){
    for(auto j = 0; j < j_max; ++j){
      for(auto k = 0; k < k_max; ++k){
        Vector3R E = Vector3R(fieldE(i,j,k,0),fieldE(i,j,k,1),fieldE(i,j,k,2));
        Vector3R J = Vector3R(fieldJ(i,j,k,0),fieldJ(i,j,k,1),fieldJ(i,j,k,2));
        potE += J.dot(E);
      }
    }
  }
  return potE;
}

Vector3R calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ, const Bounds& bounds) {
    Vector3R potE = Vector3R(0,0,0);
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
                Vector3R E = Vector3R(fieldE(i, j, k, 0), fieldE(i, j, k, 1),
                                    fieldE(i, j, k, 2));
                Vector3R J = Vector3R(fieldJ(i, j, k, 0), fieldJ(i, j, k, 1),
                                    fieldJ(i, j, k, 2));
                potE += Vector3R(J.x() * E.x(), J.y() * E.y(), J.z() * E.z() );
            }
        }
    }
    return potE;
}

// Solve Ax=b for find fieldE
void Mesh::correctE(Field3d& En, const Field3d& E, const Field3d& B,
                    Field3d& J, const double dt) {
    //zeroBoundJ(J);

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
    // TODO: use it for Field3d
    // solve_linear_system<BicgstabSolver<Field>>(
    //     A, rhs, Enew.data(), E.data());

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

void Mesh::update_Lmat(const Vector3R& coord, const Domain& domain,
                       double charge, double mass, double mpw,
                       const Field3d& fieldB, const double dt) {
    const int SMAX = SHAPE_SIZE;
    double wx, wy, wz;
    int cellLocX, cellLocY, cellLocZ, cellLocX05, cellLocY05, cellLocZ05;
    double coordLocX, coordLocY, coordLocZ;
    double coordLocX05, coordLocY05, coordLocZ05;
    int i, j, k;
    int indx, indy, indz;
    int indx05, indy05, indz05;
    alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
    alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];

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

    Vector3R B =Vector3R(0,0,0);

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

    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.5 * dt * q_m * betaI;

    const double matB[3][3] = {
        {1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
        {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
        {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    constexpr double eps = 1.e-16;

    for (int i = 0; i < SMAX; ++i) {
        for (int j = 0; j < SMAX; ++j) {
            for (int k = 0; k < SMAX; ++k) {
                // веса и индексы для (i,j,k)
                const double s1[3] = {
                    sx05[i] * sy[j] * sz[k],   // X
                    sx[i] * sy05[j] * sz[k],   // Y
                    sx[i] * sy[j] * sz05[k]    // Z
                };

                const int id1[3] = {
                    vind(cellLocX05 + i, cellLocY + j, cellLocZ + k, 0),
                    vind(cellLocX + i, cellLocY05 + j, cellLocZ + k, 1),
                    vind(cellLocX + i, cellLocY + j, cellLocZ05 + k, 2)};

                for (int i1 = 0; i1 < SMAX; ++i1) {
                    for (int j1 = 0; j1 < SMAX; ++j1) {
                        for (int k1 = 0; k1 < SMAX; ++k1) {
                            const double s2[3] = {sx05[i1] * sy[j1] * sz[k1],
                                                  sx[i1] * sy05[j1] * sz[k1],
                                                  sx[i1] * sy[j1] * sz05[k1]};

                            const int id2[3] = {
                                vind(cellLocX05 + i1, cellLocY + j1,
                                     cellLocZ + k1, 0),
                                vind(cellLocX + i1, cellLocY05 + j1,
                                     cellLocZ + k1, 1),
                                vind(cellLocX + i1, cellLocY + j1,
                                     cellLocZ05 + k1, 2)};

                            const double common = betaL;

                            // 3×3 вместо 9 копипаст
                            for (int c1 = 0; c1 < 3; ++c1) {
                                const int row = id1[c1];
                                const double w1 = s1[c1];
                                if (w1 == 0.0)
                                    continue;

                                for (int c2 = 0; c2 < 3; ++c2) {
                                    const double value =
                                        common * w1 * s2[c2] * matB[c1][c2];

                                    if (fabs(value) > eps) {
                                        LmatX[row][id2[c2]] += value;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Mesh::update_Lmat2(const Vector3R& coord, const Domain& domain,
                        double charge, double mass, double mpw,
                        const Field3d& fieldB, const double dt) {
    const int SMAX = 2;   // SHAPE_SIZE;
    alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
    alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];

    const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    const double coordLocX05 = coordLocX - 0.5;
    const double coordLocY05 = coordLocY - 0.5;
    const double coordLocZ05 = coordLocZ - 0.5;

    const int cellLocX = int(coordLocX);
    const int cellLocY = int(coordLocY);
    const int cellLocZ = int(coordLocZ);
    const int cellLocX05 = int(coordLocX05);
    const int cellLocY05 = int(coordLocY05);
    const int cellLocZ05 = int(coordLocZ05);

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

    Vector3R B = Vector3R(0.);
    // TODO: change to interpolation function
    for (int i = 0; i < SMAX; ++i) {
        const int indx = cellLocX + i;
        const int indx05 = cellLocX05 + i;
        for (int j = 0; j < SMAX; ++j) {
            const int indy = cellLocY + j;
            const int indy05 = cellLocY05 + j;
            for (int k = 0; k < SMAX; ++k) {
                const int indz = cellLocZ + k;
                const int indz05 = cellLocZ05 + k;
                const double wx = sx[i] * sy05[j] * sz05[k];
                const double wy = sx05[i] * sy[j] * sz05[k];
                const double wz = sx05[i] * sy05[j] * sz[k];
                B.x() += (wx * fieldB(indx, indy05, indz05, 0));
                B.y() += (wy * fieldB(indx05, indy, indz05, 1));
                B.z() += (wz * fieldB(indx05, indy05, indz, 2));
            }
        }
    }
    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.25 * dt * dt * q_m * betaI;

    const int blockIndex = sind(cellLocX, cellLocY, cellLocZ);
    auto& currentBlock = LmatX2[blockIndex];

    const int xOffset = cellLocX05 - cellLocX + 1;
    const int yOffset = cellLocY05 - cellLocY + 1;
    const int zOffset = cellLocZ05 - cellLocZ + 1;

    const double matB[3][3] = {
        {1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
        {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
        {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    for (int i = 0; i < SMAX; ++i) {
        for (int j = 0; j < SMAX; ++j) {
            for (int k = 0; k < SMAX; ++k) {
                const double s1[3] = {
                    sx05[i] * sy[j] * sz[k],   // X
                    sx[i] * sy05[j] * sz[k],   // Y
                    sx[i] * sy[j] * sz05[k]    // Z
                };
                const int idx1[3] = {indX(xOffset + i, j, k),
                                     indY(i, yOffset + j, k),
                                     indZ(i, j, zOffset + k)};

                for (int i1 = 0; i1 < SMAX; ++i1) {
                    for (int j1 = 0; j1 < SMAX; ++j1) {
                        for (int k1 = 0; k1 < SMAX; ++k1) {
                            const double s2[3] = {sx05[i1] * sy[j1] * sz[k1],
                                                  sx[i1] * sy05[j1] * sz[k1],
                                                  sx[i1] * sy[j1] * sz05[k1]};
                            const int idx2[3] = {indX(xOffset + i1, j1, k1),
                                                 indY(i1, yOffset + j1, k1),
                                                 indZ(i1, j1, zOffset + k1)};

                            for (int c1 = 0; c1 < 3; ++c1) {
                                const int rowIndex = idx1[c1];
                                for (int c2 = 0; c2 < 3; ++c2) {
                                    const int colIndex = idx2[c2];
                                        currentBlock(rowIndex, colIndex,
                                                     c1 * 3 + c2) +=
                                        betaL * s1[c1] * s2[c2] * matB[c1][c2];
                                }
                            }
                        }   // k1
                    }   // j1
                }   // i1
            }   // k
        }   // j
    }   // i
}

void Mesh::update_LmatNGP(const Vector3R& coord, const Domain& domain,
                          double charge, double mass, double mpw,
                          const Field3d& fieldB, const double dt) {
    const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    const double coordLocX05 = coordLocX - 0.5;
    const double coordLocY05 = coordLocY - 0.5;
    const double coordLocZ05 = coordLocZ - 0.5;

    const int cellLocX = ngp(coordLocX);
    const int cellLocY = ngp(coordLocY);
    const int cellLocZ = ngp(coordLocZ);
    const int cellLocX05 = ngp(coordLocX05);
    const int cellLocY05 = ngp(coordLocY05);
    const int cellLocZ05 = ngp(coordLocZ05);

    const int indx = vind(cellLocX05, cellLocY, cellLocZ, 0);
    const int indy = vind(cellLocX, cellLocY05, cellLocZ, 1);
    const int indz = vind(cellLocX, cellLocY, cellLocZ05, 2);
    Vector3R B = Vector3R(0.);

    B.x() = fieldB(cellLocX, cellLocY05, cellLocZ05, 0);
    B.y() = fieldB(cellLocX05, cellLocY, cellLocZ05, 1);
    B.z() = fieldB(cellLocX05, cellLocY05, cellLocZ, 2);

    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.5 * dt * q_m * betaI;

    const double matB[3][3] = {
        {1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
        {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
        {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    constexpr double eps = 1.e-16;
    const double common = betaL;

    const int id[3] = {indx, indy, indz};

    for (int c1 = 0; c1 < 3; ++c1) {
        const int row = id[c1];
        for (int c2 = 0; c2 < 3; ++c2) {
            const double value = common * matB[c1][c2];
            if (fabs(value) > eps) {
                LmatX[row][id[c2]] += value;
            }
        }
    }
}

void Mesh::update_Lmat2_NGP(const Vector3R& coord, const Domain& domain,
                            double charge, double mass, double mpw,
                            const Field3d& fieldB, const double dt) {
    Vector3R B = Vector3R(0.);
    const double coordLocX = coord.x() / domain.cell_size().x() + GHOST_CELLS;
    const double coordLocY = coord.y() / domain.cell_size().y() + GHOST_CELLS;
    const double coordLocZ = coord.z() / domain.cell_size().z() + GHOST_CELLS;
    const double coordLocX05 = coordLocX - 0.5;
    const double coordLocY05 = coordLocY - 0.5;
    const double coordLocZ05 = coordLocZ - 0.5;

    const int cellLocX = ngp(coordLocX);
    const int cellLocY = ngp(coordLocY);
    const int cellLocZ = ngp(coordLocZ);
    const int cellLocX05 = ngp(coordLocX05);
    const int cellLocY05 = ngp(coordLocY05);
    const int cellLocZ05 = ngp(coordLocZ05);

    B.x() = fieldB(cellLocX, cellLocY05, cellLocZ05, 0);
    B.y() = fieldB(cellLocX05, cellLocY, cellLocZ05, 1);
    B.z() = fieldB(cellLocX05, cellLocY05, cellLocZ, 2);

    const double q_m = charge / mass;
    const Vector3R b = 0.5 * dt * q_m * B;

    const double betaI = mpw * charge / (1.0 + b.squared());
    const double betaL = 0.5 * dt * q_m * betaI;

    const int blockIndex = sind(cellLocX, cellLocY, cellLocZ);
    auto& currentBlock = LmatX2[blockIndex];

    const int xOffset = cellLocX05 - cellLocX + 1;
    const int yOffset = cellLocY05 - cellLocY + 1;
    const int zOffset = cellLocZ05 - cellLocZ + 1;

    const int indx = BlockDimsNGP::indX(0, yOffset, zOffset);
    const int indy = BlockDimsNGP::indY(xOffset, 0, zOffset);
    const int indz = BlockDimsNGP::indZ(xOffset, yOffset, 0);

    const double matB[3][3] = {
        {1.0 + b.x() * b.x(), +b.z() + b.x() * b.y(), -b.y() + b.x() * b.z()},
        {-b.z() + b.y() * b.x(), 1.0 + b.y() * b.y(), +b.x() + b.y() * b.z()},
        {+b.y() + b.z() * b.x(), -b.x() + b.z() * b.y(), 1.0 + b.z() * b.z()}};

    const double common = betaL;
    const int id[3] = {indx, indy, indz};

    for (int c1 = 0; c1 < 3; ++c1) {
        for (int c2 = 0; c2 < 3; ++c2) {
            currentBlock(id[c1], id[c2], c1 * 3 + c2) += common * matB[c1][c2];
        }
    }
}

void Mesh::apply_periodic_boundaries(std::vector<IndexMap>& LmatX) {
    const auto size = Vector3I(xSize, ySize, zSize);
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

    const auto size = Vector3I(xSize, ySize, zSize);
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
            for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
                auto ind2 = it->first;
                setValuesZero(it->second, domain, i, ind2);
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

    const auto size = Vector3I(xSize, ySize, zSize);
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


void Mesh::apply_open_boundaries(Operator& LmatX, Domain& domain) {

    const auto size = Vector3I(xSize, ySize, zSize);
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
