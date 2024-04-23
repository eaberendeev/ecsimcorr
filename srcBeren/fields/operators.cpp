#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "bounds.h"

// matrix ColMajor
// (0,0) (0,1) (0,2)
// (1,0) (1,1) (1,2)
// (2,0) (2,1) (2,2)

// Ex(i+/2,j,k), Ey(i,j+1/2,k), Ez(i,j,k+1/2) 
// Bx(i,j+1/2,k+1/2), By(i+1/2,j,k+1/2), Bz(i+/2,j+1/2,k)

 void Mesh::stencil_Imat(const Domain &domain)
 { 
 // !!!!! needs bound condition and if cases!!!!!!
  std::vector<Trip> trips;
  const auto size = domain.size();
  int totalSize = size.x()*size.y()*size.z();
  trips.reserve(totalSize);
    
  for(int i = 0; i < size.x(); i++){
    for(int j = 0; j < size.y(); j++){
      for(int k = 0; k < size.z(); k++){

        // i,j,k
        trips.push_back(Trip(vind(i,j,k,0),vind(i,j,k,0),1.0));

        // i,j,k
        trips.push_back(Trip(vind(i,j,k,1),vind(i,j,k,1),1.0));
          
        // i,j,k
        trips.push_back(Trip(vind(i,j,k,2),vind(i,j,k,2),1.0));
          
      }
    }
  }
  std::cout << size << " " << 3*totalSize << " " << Imat.rows() << " " << Imat.cols() <<  "\n";
  Imat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_Lmat(const Domain &domain) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    const int totalSize = 20*size.x()*size.y()*size.z()*SHAPE_SIZE*SHAPE_SIZE*SHAPE_SIZE;
    std::cout << totalSize << "\n";
    trips.reserve(totalSize);
    int ind2;
    double value;
    
    for(int i = 0; i < 3*(size.x()*size.y()*size.z() ) ; i++){
      for (auto it=LmatX[i].begin(); it!=LmatX[i].end(); ++it){
        ind2 = it->first;
        value = it->second;
	if(fabs(value) > 1.e-16)
          trips.push_back(Trip(i,ind2, value ));
      }
    }
    
    Lmat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlB(const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x()*size.y()*size.z()*12;
    trips.reserve(totalSize);
    double val;
    int vindx, vindy, vindz;
    int im,jm,km;
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();

    for(int i = 0; i < size.x(); i++){
      for(int j = 0; j < size.y(); j++){
        for(int k = 0; k < size.z(); k++){

          if( i!= 0){
            im = i - 1;
          } else{
			      if(isPeriodicX){
              im = size.x() - 4;
            } else {
              continue;
            }
          }
          if( j!= 0){
            jm = j - 1;
          } else{
			      if(isPeriodicY){
              jm = size.y() - 4;
            } else {
              continue;
            }
          }
          if( k!= 0){
            km = k - 1;
          } else{
			      if(isPeriodicZ){
              km = size.z() - 4;
            } else {
              continue;
            }
          }

          vindx = vind(i,j,k,0);
          vindy = vind(i,j,k,1);
          vindz = vind(i,j,k,2);

          // (x)[i+1/2,j,k] 
          // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
          val = 1.0 / dy;
          trips.push_back(Trip(vindx,vind(i , j , k , 2), val));
          trips.push_back(Trip(vindx,vind(i , jm, k , 2),-val));
          // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
          val = -1.0 / dz;
          trips.push_back(Trip(vindx,vind(i , j , k , 1), val));
          trips.push_back(Trip(vindx,vind(i , j , km, 1),-val));

          // (y)[i,j+1/2,k] 
          // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
          val = 1.0 / dz;
          trips.push_back(Trip(vindy,vind(i , j , k , 0), val));
          trips.push_back(Trip(vindy,vind(i , j , km, 0),-val));
          // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
          val = -1.0 / dx;
          trips.push_back(Trip(vindy,vind(i , j , k, 2), val));
          trips.push_back(Trip(vindy,vind(im, j , k, 2),-val));

          // (z)[i,j,k+1/2] 
          // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
          val = 1.0 / dx;
          trips.push_back(Trip(vindz,vind(i , j , k , 1), val));
          trips.push_back(Trip(vindz,vind(im, j , k , 1),-val));
          // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
          val = -1.0 / dy;
          trips.push_back(Trip(vindz,vind(i , j , k , 0), val));
          trips.push_back(Trip(vindz,vind(i , jm, k , 0),-val));

        }
      }
    }
    curlB.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlE(const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = fieldE.size();
    int totalSize = size.x()*size.y()*size.z()*12;
    trips.reserve(totalSize);
    double val;
    int vindx, vindy, vindz;
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();
    int ip,jp,kp;

    for(int i = 0; i < size.x(); i++){
      for(int j = 0; j < size.y(); j++){
        for(int k = 0; k < size.z(); k++){

          if( i!= size.x() - 1){
            ip = i + 1;
          } else{
			      if(isPeriodicX){
              ip = 3;
            } else {
              continue;
            }
          }
          if( j!= size.y() - 1){
            jp = j + 1;
          } else{
			      if(isPeriodicY){
              jp = 3;
            } else {
              continue;
            }
          }
          if( k!= size.z() - 1){
            kp = k + 1;
          } else{
			      if(isPeriodicZ){
              kp = 3;
            } else {
              continue;
            }
          }

          vindx = vind(i,j,k,0);
          vindy = vind(i,j,k,1);
          vindz = vind(i,j,k,2);

          // (x)[i,j+1/2,k+1/2] 
          // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
          val = 1.0 / dy;
          trips.push_back(Trip(vindx,vind(i , jp, k , 2), val));
          trips.push_back(Trip(vindx,vind(i , j,  k , 2),-val));
          // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
          val = -1.0 / dz;
          trips.push_back(Trip(vindx,vind(i , j , kp, 1), val));
          trips.push_back(Trip(vindx,vind(i , j , k , 1),-val));

          // (y)[i+1/2,j,k+1/2] 
          // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
          val = 1.0 / dz;
          trips.push_back(Trip(vindy,vind(i , j ,kp, 0), val));
          trips.push_back(Trip(vindy,vind(i , j ,k , 0),-val));
          // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
          val = -1.0 / dx;
          trips.push_back(Trip(vindy,vind(ip, j , k , 2), val));
          trips.push_back(Trip(vindy,vind(i , j , k , 2),-val));

          // (z)[i+1/2,j+1/2,k] 
          // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
          val = 1.0 / dx;
          trips.push_back(Trip(vindz,vind(ip, j , k , 1), val));
          trips.push_back(Trip(vindz,vind(i , j , k , 1),-val));
          // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
          val = -1.0 / dy;
          trips.push_back(Trip(vindz,vind(i , jp, k , 0), val));
          trips.push_back(Trip(vindz,vind(i , j , k , 0),-val));

        }
      }
    }
    curlE.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_divE(const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = fieldE.size();
    int totalSize = size.x()*size.y()*size.z()*6;
    trips.reserve(totalSize);
    double val;
    int sindx;
    int im,jm,km;
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();
    for(int i = 0; i < size.x(); i++){
      for(int j = 0; j < size.y(); j++){
        for(int k = 0; k < size.z(); k++){


          if( i!= 0){
            im = i - 1;
          } else{
			      if(isPeriodicX){
              im = size.x() - 4;
            } else {
              continue;
            }
          }
          if( j!= 0){
            jm = j - 1;
          } else{
			      if(isPeriodicY){
              jm = size.y() - 4;
            } else {
              continue;
            }
          }
          if( k!= 0){
            km = k - 1;
          } else{
			      if(isPeriodicZ){
              km = size.z() - 4;
            } else {
              continue;
            }
          }
          
          sindx = sind(i,j,k);

          // [i,j,k] 
          // ( Ex[i+1/2,j,k] - Ex[i-1,j,k] ) / dx
          val = 1.0 / dx;
          trips.push_back(Trip(sindx,vind(i , j, k , 0), val));
          trips.push_back(Trip(sindx,vind(im, j, k , 0),-val));
          // ( Ex[i,j+1/2,k] - Ex[i,j-1/2,k] ) / dy
          val = 1.0 / dy;
          trips.push_back(Trip(sindx,vind(i , j , k , 1), val));
          trips.push_back(Trip(sindx,vind(i , jm, k , 1),-val));
          // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
          val = 1.0 / dz;
          trips.push_back(Trip(sindx,vind(i , j , k , 2), val));
          trips.push_back(Trip(sindx,vind(i , j , km, 2),-val));

        }
      }
    }
    divE.setFromTriplets(trips.begin(), trips.end());
}
