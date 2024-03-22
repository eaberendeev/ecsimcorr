#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "solverSLE.h"
#include "bounds.h"

Mesh::Mesh(const World& world) : 
          Lmat(world.region.total_size()*3,world.region.total_size()*3 ),
          Lmat2(world.region.total_size()*3,world.region.total_size()*3 ),
          Mmat(world.region.total_size()*3,world.region.total_size()*3 ),
          Imat(world.region.total_size()*3,world.region.total_size()*3 ),
          curlE(world.region.total_size()*3,world.region.total_size()*3 ),
          curlB(world.region.total_size()*3,world.region.total_size()*3 ),
          fieldE(world.region.numNodes,3),
          fieldEn(world.region.numNodes,3),
          fieldEp(world.region.numNodes,3),
          fieldB(world.region.numNodes,3),
          fieldJp(world.region.numNodes,3),
          fieldJp_full(world.region.numNodes,3),
          fieldJe(world.region.numNodes,3),
          fieldB0(world.region.numNodes,3), 
          fieldBInit(world.region.numNodes,3),
          LmatX(world.region.total_size()*3), 
          chargeDensityOld(world.region.numNodes,1), 
          chargeDensity(world.region.numNodes,1), 
          divE(world.region.total_size(),world.region.total_size()*3), 
          _world(world)
{
    
       
  
  set_fields();
  
  if(RECOVERY > 0) {
    read_fields_from_recovery();
  }
    _size1 = world.region.numNodes.x();
    _size2 = world.region.numNodes.y();
    _size3 = world.region.numNodes.z();

    stencil_Imat();
    stencil_curlE();
    stencil_curlB();
    stencil_divE();

    Mmat = -0.25*Dt*Dt*curlB*curlE;

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


void Mesh::write_field_to_file(const std::string& dataName, Field3d& field){
    std::ofstream file_bin(dataName, std::ios::out | std::ios::binary);
    int size_x, size_y, size_z;
    int3 size = field.size();
    size_x = size.x();
    size_y = size.y();
    size_z = size.z();
    
    file_bin.write((char*) &size_x, sizeof(size_x));
    file_bin.write((char*) &size_y, sizeof(size_y));
    file_bin.write((char*) &size_z, sizeof(size_z));

  int dataSize = 3*size_x*size_y*size_z;
  std::vector<double> data3d(dataSize);
  for(auto i = 0; i < dataSize; i++ ){
            data3d[i] = field(i);
  }

  file_bin.write((char*) &data3d[0], dataSize*sizeof(data3d[0]));

  file_bin.close();

}

void Mesh::read_field_from_file(const std::string& dataName, Field3d& field){
    std::ifstream file_bin(dataName, std::ios::in | std::ios::binary);
    int size_x, size_y, size_z;
    file_bin.read((char*) &size_x, sizeof(size_x));
    file_bin.read((char*) &size_y, sizeof(size_y));
    file_bin.read((char*) &size_z, sizeof(size_z));
    std::cout<< size_x << " " << size_y << " "<< size_z << "\n";
    int3 size = field.size();
  if(size_x != size.x() || size_y != size.y() || size_z != size.z()) {
    std::cout << "Invalid reading from recovery! Field size is invalid! \n";
    exit(0);
  }

  int dataSize = 3*size_x*size_y*size_z;
  std::vector<double> data3d(dataSize);
  file_bin.read((char*) &data3d[0], dataSize*sizeof(data3d[0]));
  for(auto i = 0; i < dataSize; i++ ){
            field(i) = data3d[i];

}
  file_bin.close();

}

void Mesh::read_fields_from_recovery(){
    std::string fname = "..//Recovery//Fields//FieldE.backup";
    read_field_from_file(fname, fieldEn);
    fname = "..//Recovery//Fields//FieldB.backup";
    read_field_from_file(fname, fieldB);
} 

void Mesh::write_fields_to_recovery(int timestep){
    if( RecTimeStep < 0){
        return;
    }
    if( timestep % RecTimeStep != 0){
        return;
    }

    std::cout << "Backup mesh in " << timestep << "\n";
    std::string fname = ".//Recovery//Fields//FieldE.backup";
    write_field_to_file(fname, fieldEn);
    fname = ".//Recovery//Fields//FieldB.backup";
    write_field_to_file(fname, fieldB);
} 

void Mesh::set_fields(){
	set_uniform_fields();

  fieldEn = fieldE;
  fieldB = fieldBInit;
  fieldB0 = fieldB;

  fieldJp = 0.;
  fieldJe = 0.;
  fieldEp = 0.;
} 


void Mesh::set_uniform_fields(){
  auto size = fieldE.size();
  fieldE = 0.;
  for( auto i = 0; i < size.x(); ++i){
    for( auto j = 0; j < size.y(); ++j){
      for( auto k = 0; k < size.z(); ++k){
        fieldBInit(i,j,k,0) = BUniform[0];
        fieldBInit(i,j,k,1) = BUniform[1];
        fieldBInit(i,j,k,2) = BUniform[2];
      }
    }
  }
  
}

void Mesh::prepare()
{
  fieldE = fieldEn;
  fieldB0 = fieldB;
  fieldJp = 0.;
  fieldJe = 0.;

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
  if(isPeriodicX){
    i_max -= 3;
  }
  if(isPeriodicY){
    j_max -= 3;
  }
  if(isPeriodicZ){
    k_max -= 3;
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

double Mesh::calc_JE(const Field3d& fieldE,const Field3d& fieldJ) const{
  double potE = 0;
  int i_max = fieldE.size().x();
  int j_max = fieldE.size().y();
  int k_max = fieldE.size().z();
  
  if(isPeriodicX){
    i_max -= 3;
  }
  if(isPeriodicY){
    j_max -= 3;
  }
  if(isPeriodicZ){
    k_max -= 3;
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

double3 Mesh::calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ) const {
    double3 potE = double3(0,0,0);
    int i_max = fieldE.size().x();
    int j_max = fieldE.size().y();
    int k_max = fieldE.size().z();

    if (isPeriodicX) {
        i_max -= 3;
    }
    if (isPeriodicY) {
        j_max -= 3;
    }
    if (isPeriodicZ) {
        k_max -= 3;
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


void Mesh::correctE()
{
  fieldB.data() -= fieldBInit.data();

	  static Field rhs;
    rhs = fieldE.data() - Dt*fieldJe.data() + Dt*curlB*fieldB.data() + Mmat*fieldE.data();
	  static Operator A = Imat - Mmat;
    solve_SLE(A, rhs, fieldEn.data(), fieldE.data());
    std::cout<< "Solver2 error = "<< (A*fieldEn.data() - rhs).norm() << "\n";

  fieldB.data() += fieldBInit.data();

}

void Mesh::predictE() {
    fieldB.data() -= fieldBInit.data();

    Lmat2 = Lmat;

    Lmat = Mmat - Lmat;

    static Field rhs;
    rhs = fieldE.data() - Dt * fieldJp.data() + Dt * curlB * fieldB.data() +
          Lmat * fieldE.data();
    Operator A = Imat - Lmat;
    solve_SLE(A, rhs, fieldEp.data(), fieldE.data());

    std::cout << "Solver1 error = " << (A * fieldEp.data() - rhs).norm()
              << "\n";

    fieldB.data() += fieldBInit.data();
}

void Mesh::fdtd_explicit()
{
  fieldB.data() -= fieldBInit.data();

  fieldE.data() += 0.5*Dt*curlB*fieldB.data() - 0.5*Dt*fieldJp.data();
  fieldB.data() -= 0.5*Dt*curlE*fieldE.data();
  fieldB.data() += fieldBInit.data();

}

void Mesh::computeB()
{
    fieldB.data() -= 0.5*Dt*curlE*(fieldE.data() + fieldEn.data() );
}

void Mesh::update_Lmat( const double3& coord, double charge, double mass, double mpw)
{
    const int SMAX = SHAPE_SIZE;
    double3 B,h;
    double alpha;
    double Ap ;
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

        coordLocX = coord.x() / Dx + CELLS_SHIFT;
        coordLocY = coord.y() / Dy + CELLS_SHIFT;
        coordLocZ = coord.z() / Dz + CELLS_SHIFT;
        coordLocX05 = coord.x() / Dx + CELLS_SHIFT - 0.5;
        coordLocY05 = coord.y() / Dy + CELLS_SHIFT - 0.5;
        coordLocZ05 = coord.z() / Dz + CELLS_SHIFT - 0.5;

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
        h = unit(B);

        alpha = 0.5*Dt*charge*mag(B) / mass;
        Ap = 0.25*Dt*Dt*mpw*charge*charge / mass / (1+alpha*alpha); 
        
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

void Mesh::glue_Lmat_bound()
{
    const auto size = fieldE.size();
    const int overlap = 3;
    const int last_indx = size.x() - overlap;
    const int last_indy = size.y() - overlap;
    const int last_indz = size.z() - overlap;

  if(isPeriodicX){
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( ix < overlap ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(last_indx + ix,iy,iz,id);
            
            if(ix1 < overlap){
              auto indBound2 = vind(last_indx + ix1,iy1,iz1,id1);
              LmatX[indBound][indBound2] += value;
            } else {
              LmatX[indBound][ind2] += value;              
            }
          }
        } 
    }
  }

  if(isPeriodicY){
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( iy < overlap ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(ix,last_indy + iy,iz,id);
            if(iy1 < overlap){
              auto indBound2 = vind(ix1,last_indy + iy1,iz1,id1);
              LmatX[indBound][indBound2] += value;
            } else {
              LmatX[indBound][ind2] += value;              
            }
          }
        }
    }
  }

  if(isPeriodicZ){
#pragma omp parallel for
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
        auto ix = pos_vind(i,0); 
        auto iy = pos_vind(i,1);
        auto iz = pos_vind(i,2); 
        auto id = pos_vind(i,3); 
        
        if( iz < overlap ){
          for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
            auto ind2 = it->first;
            auto value = it->second;
            auto ix1 = pos_vind(ind2,0); 
            auto iy1 = pos_vind(ind2,1);
            auto iz1 = pos_vind(ind2,2); 
            auto id1 = pos_vind(ind2,3); 
            auto indBound = vind(ix,iy,last_indz + iz,id);
            if(iz1 < overlap){
              auto indBound2 = vind(ix1,iy1,last_indz + iz1,id1);
              LmatX[indBound][indBound2] += value;
            } else {
              LmatX[indBound][ind2] += value;              
            }
          }
        }
    }
  }
  
  if(isPeriodicX){
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

  if(isPeriodicY){
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

  if(isPeriodicZ){
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

void Mesh::make_periodic_border_with_add(Field3d &field ){

  auto size = field.size();

  double overlap = 3;
  auto i_max = size.x() - overlap; 
  auto j_max = size.y() - overlap; 
  auto k_max = size.z() - overlap; 

if(isPeriodicX){
  for(auto i = 0; i < overlap; ++i){
    for(auto j = 0; j < size.y(); ++j){
      for(auto k = 0; k < size.z(); ++k){
        for(auto dim = 0; dim < 3; dim++){
          field(i,j,k,dim) += field(i + i_max,j,k,dim);
          field(i + i_max,j,k,dim) = field(i,j,k,dim);
        } 
      }
    }
  }
}

if(isPeriodicY){
  for(auto i = 0; i < size.x(); ++i){
    for(auto j = 0; j < overlap; ++j){
      for(auto k = 0; k < size.z(); ++k){
        for(auto dim = 0; dim < 3; dim++){
          field(i,j,k,dim) += field(i,j + j_max,k,dim);
          field(i,j + j_max,k,dim) = field(i,j,k,dim);
        } 
      }
    }
  }
  }
if(isPeriodicZ){
  for(auto i = 0; i < size.x(); ++i){
    for(auto j = 0; j < size.y(); ++j){
      for(auto k = 0; k < overlap; ++k){
        for(auto dim = 0; dim < 3; dim++){
          field(i,j,k,dim) += field(i,j,k + k_max,dim);
          field(i,j,k + k_max,dim) = field(i,j,k,dim);
        } 
      }
    }
  }
}
}

double3 Mesh::get_fieldE_in_cell(int i, int j,int k)  const{
  double3 E;
  E.x() = 0.25 * (fieldE(i,j,  k,0) + fieldE(i,j  ,k+1,0) 
                + fieldE(i,j+1,k,0) + fieldE(i,j+1,k+1,0) );
  
  E.y() = 0.25 * (fieldE(i,  j,k,1) + fieldE(i,  j,k+1,1) + 
                  fieldE(i+1,j,k,1) + fieldE(i+1,j,k+1,1) );
  
  E.z() = 0.125 * (fieldE(i,  j,  k,2) + fieldE(i,  j,  k+1,2) + 
                   fieldE(i,  j+1,k,2) + fieldE(i,  j+1,k+1,2) +
                   fieldE(i+1,j,  k,2) + fieldE(i+1,j,  k+1,2) + 
                   fieldE(i+1,j+1,k,2) + fieldE(i+1,j+1,k+1,2) );
  return E;
}

double3 Mesh::get_fieldB_in_cell(int i, int j,int k)  const{
  double3 B;
  B.x() = 0.25 * (fieldB(i,  j,k,0) + fieldB(i,  j,k+1,0) +
                  fieldB(i+1,j,k,0) + fieldB(i+1,j,k+1,0) );
  
  B.y() = 0.25 * (fieldB(i,j,  k,1) + fieldB(i,j,  k+1,1) +
                  fieldB(i,j+1,k,1) + fieldB(i,j+1,k+1,1) );
  
  B.z() = fieldB(i,j,k,2);
  return B;
}


double3 get_fieldE_in_pos(const Field3d& fieldE,const double3& r) {

    double3 E;
    int indx, indy,indz,indx1, indy1,indz1;

    double xx,yy,zz;

    double sx0,sy0,sz0,sdx0,sdy0,sdz0;
    double sx1,sy1,sz1,sdx1,sdy1,sdz1;

    const double rdx = 1. / Dx;
    const double rdy = 1. / Dy;
    const double rdz = 1. / Dz;    
    xx = r.x() * rdx;
    yy = r.y() * rdy;
    zz = r.z() * rdz;
    
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
double3 get_fieldB_in_pos(const Field3d& fieldB, const double3& r) {

    double3 B;
    int indx, indy,indz,indx1, indy1,indz1;

    double xx,yy,zz;

    double sx0,sy0,sz0,sdx0,sdy0,sdz0;
    double sx1,sy1,sz1,sdx1,sdy1,sdz1;

    const double rdx = 1. / Dx;
    const double rdy = 1. / Dy;
    const double rdz = 1. / Dz;    
    xx = r.x() * rdx;
    yy = r.y() * rdy;
    zz = r.z() * rdz;
    
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

void get_fields_in_pos(const Field3d& fieldE,const Field3d& fieldB, const double3& r, double3& locE, double3 &locB) {

    int indx, indy,indz,indx1, indy1,indz1;

    double xx,yy,zz;

    double sx0,sy0,sz0,sdx0,sdy0,sdz0;
    double sx1,sy1,sz1,sdx1,sdy1,sdz1;

    const double rdx = 1. / Dx;
    const double rdy = 1. / Dy;
    const double rdz = 1. / Dz;

    xx = r.x() * rdx;
    yy = r.y() * rdy;
    zz = r.z() * rdz;
    
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
    locE.x() = sdx0 * ( sy0 * ( sz0 * fieldE(indx1,indy,indz,0) + sz1 * fieldE(indx1,indy,indz+1,0) ) 
                + sy1 * ( sz0 * fieldE(indx1,indy+1,indz,0) + sz1 * fieldE(indx1,indy+1,indz+1,0) ) ) 
        + sdx1 * ( sy0 * ( sz0 * fieldE(indx1+1,indy,indz,0) + sz1 * fieldE(indx1+1,indy,indz+1,0) ) 
                + sy1 * ( sz0 * fieldE(indx1+1,indy+1,indz,0) + sz1 * fieldE(indx1+1,indy+1,indz+1,0) ) );

    locE.y() = sx0 * ( sdy0 * ( sz0 * fieldE(indx,indy1,indz,1) + sz1 * fieldE(indx,indy1,indz+1,1) ) 
                + sdy1 * ( sz0 * fieldE(indx,indy1+1,indz,1) + sz1 * fieldE(indx,indy1+1,indz+1,1) ) ) 
        + sx1 * ( sdy0 * ( sz0 * fieldE(indx+1,indy1,indz,1) + sz1 * fieldE(indx+1,indy1,indz+1,1) ) 
                + sdy1 * ( sz0 * fieldE(indx+1,indy1+1,indz,1) + sz1 * fieldE(indx+1,indy1+1,indz+1,1) ) );

    locE.z() = sx0 * ( sy0 * ( sdz0 * fieldE(indx,indy,indz1,2) + sdz1 * fieldE(indx,indy,indz1+1,2) ) 
                + sy1 * ( sdz0 * fieldE(indx,indy+1,indz1,2) + sdz1 * fieldE(indx,indy+1,indz1+1,2) ) ) 
        + sx1 * ( sy0 * ( sdz0 * fieldE(indx+1,indy,indz1,2) + sdz1 * fieldE(indx+1,indy,indz1+1,2) ) 
                + sy1 * ( sdz0 * fieldE(indx+1,indy+1,indz1,2) + sdz1 * fieldE(indx+1,indy+1,indz1+1,2) ) );

  locB.x() = sx0 * ( sdy0 * ( sdz0 * fieldB(indx,indy1,indz1,0) + sdz1 * fieldB(indx,indy1,indz1+1,0) ) 
               + sdy1 * ( sdz0 * fieldB(indx,indy1+1,indz1,0) + sdz1 * fieldB(indx,indy1+1,indz1+1,0) ) ) 
       + sx1 * ( sdy0 * ( sdz0 * fieldB(indx+1,indy1,indz1,0) + sdz1 * fieldB(indx+1,indy1,indz1+1,0) ) 
               + sdy1 * ( sdz0 * fieldB(indx+1,indy1+1,indz1,0) + sdz1 * fieldB(indx+1,indy1+1,indz1+1,0) ) );

  locB.y() = sdx0 * ( sy0 * ( sdz0 * fieldB(indx1,indy,indz1,1) + sdz1 * fieldB(indx1,indy,indz1+1,1) ) 
              + sy1 * ( sdz0 * fieldB(indx1,indy+1,indz1,1) + sdz1 * fieldB(indx1,indy+1,indz1+1,1) ) ) 
      + sdx1 * ( sy0 * ( sdz0 * fieldB(indx1+1,indy,indz1,1) + sdz1 * fieldB(indx1+1,indy,indz1+1,1) ) 
              + sy1 * ( sdz0 * fieldB(indx1+1,indy+1,indz1,1) + sdz1 * fieldB(indx1+1,indy+1,indz1+1,1) ) );

  locB.z() = sdx0 * ( sdy0 * ( sz0 * fieldB(indx1,indy1,indz,2) + sz1 * fieldB(indx1,indy1,indz+1,2) ) 
              + sdy1 * ( sz0 * fieldB(indx1,indy1+1,indz,2) + sz1 * fieldB(indx1,indy1+1,indz+1,2) ) ) 
      + sdx1 * ( sdy0 * ( sz0 * fieldB(indx1+1,indy1,indz,2) + sz1 * fieldB(indx1+1,indy1,indz+1,2) ) 
              + sdy1 * ( sz0 * fieldB(indx1+1,indy1+1,indz,2) + sz1 * fieldB(indx1+1,indy1+1,indz+1,2) ) );

}
