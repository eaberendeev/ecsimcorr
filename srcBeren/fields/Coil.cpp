#include "Coil.h"

double CoilsArray::get_integ_r(double z,double r,double R){
    double sum = 0;
    double znam;
    for (auto i = 0; i < N; i++){
        znam =  R*R + z*z + r*r - 2*R*r*cs[i];
        sum += hp*(cs[i] / (znam*sqrt(znam) ) );
    }
    return sum;
}
double CoilsArray::get_integ_z(double z,double r,double R){
    double sum = 0;
    double znam;

    for (auto i = 0; i < N; i++){
        znam =  R*R + z*z + r*r - 2*R*r*cs[i];
        sum += hp*((R- r*cs[i] )/ (znam*sqrt(znam) ) );
    }
    return sum;
}



double CoilsArray::get_Bz(double z,double r){
    double Bz = 0;
    for (const auto& coil : coils){
    	double zc = z - coil.z0;
    	Bz += coil.I*coil.R*get_integ_z(zc,r,coil.R);
	}
    return Bz;
}
double CoilsArray::get_Br(double z,double r){
    double Br = 0;
    for (const auto& coil : coils){
      double zc = z - coil.z0;
      Br += coil.I*coil.R*zc*get_integ_r(zc,r,coil.R);
  }
    return Br;
}


void set_coils( Field3d& fieldB,const World& world){
  auto size_x = fieldB.size().x();
  auto size_y = fieldB.size().y();
  auto size_z = fieldB.size().z();
  double center_y = 0.5*world.regionGlob.numCells.y()*Dy;
  double center_z = 0.5*world.regionGlob.numCells.z()*Dz;
  double yy,rr,zz;
  double Bry,Brz,Bz;
  CoilsArray coils;



        //std::cout<< coil.z0 << " " <<coil.R << " "<<amp  <<"\n";
    #pragma omp parallel for collapse(1) private(yy,rr,zz,Bry,Brz,Bz)
    for( auto i = 0; i < size_x; ++i){
        double xx = i*Dx+world.region.origin - Dx*CELLS_SHIFT;
      for( auto j = 0; j < size_y; ++j){
        for( auto k = 0; k < size_z; ++k){
      	    yy = (j+0.5)*Dy - center_y - Dy*CELLS_SHIFT;
            zz = (k+0.5)*Dz - center_z - Dz*CELLS_SHIFT; 
            rr = sqrt(zz*zz+yy*yy);
            Bz = coils.get_Bz(xx,rr);
            fieldB(i,j,k,0) += Bz;

            yy = j*Dy - center_y- Dy*CELLS_SHIFT;
            zz = (k+0.5)*Dz - center_z- Dz*CELLS_SHIFT;
            rr = sqrt(zz*zz+yy*yy);            
            Bry = coils.get_Br(xx+0.5*Dx,rr);
            fieldB(i,j,k,1) += Bry*yy/rr;            
            
            yy = (j+0.5)*Dy - center_y- Dy*CELLS_SHIFT;
            zz = k*Dz - center_z- Dz*CELLS_SHIFT; 
            rr = sqrt(zz*zz+yy*yy);     
            Brz = coils.get_Br(xx+0.5*Dx,rr);          
            fieldB(i,j,k,2) += Brz*zz/rr;
        }
      }
    }

    std::cout<<"Coils configuration set successfull!\n";
}
