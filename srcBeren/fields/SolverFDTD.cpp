#include "World.h"
#include "Mesh.h"
#include "SolverFDTD.h"
#include "Vec.h"
//static void update_fieldsB(const Field3d& fieldE, Field3d& fieldB){
    // int i, j, k;
    // const double dtp = 0.5 * Dt;
    // const double rdx = 1. / Dx;
    // const double rdy = 1. / Dy;
    // const double rdz = 1. / Dz;
    // const int size_d1 = fieldB.size().x();
    // const int size_d2 = fieldB.size().y();
    // const int size_d3 = fieldB.size().z();
	
/*
    for( i = 0; i < size_d1 - 1; ++i ){
        for( j = 0; j < size_d2 - 1; ++j ){
            for( k = 0; k < size_d3 - 1; ++k ){
    			fieldB(i,j,k).x() -= dtp * (  rdy * (fieldE(i,j+1,k).z() - fieldE(i,j,k).z() ) 
                                              - rdz * (fieldE(i,j,k+1).y() - fieldE(i,j,k).y() )  );


    			fieldB(i,j,k).y() -= dtp * ( rdz * (fieldE(i,j,k+1).x() - fieldE(i,j,k).x() )
                                           - rdx * (fieldE(i+1,j,k).z() - fieldE(i,j,k).z() ) );


    			fieldB(i,j,k).z() -= dtp * ( rdx * (fieldE(i+1,j,k).y() - fieldE(i,j,k).y() )
                                           - rdy * (fieldE(i,j+1,k).x() - fieldE(i,j,k).x() ) );
            }
		}
    }

    i = size_d1 - 1;
    for( j = 0; j < size_d2 - 1; ++j){
        for( k = 0; k < size_d3 - 1; ++k){

	       	fieldB(i,j,k).x() -= dtp * (  rdy * (fieldE(i,j+1,k).z() - fieldE(i,j,k).z() ) 
                                              - rdz * (fieldE(i,j,k+1).y() - fieldE(i,j,k).y() )  );
        }
    }

	j = size_d2 - 1;
    for(i = 0; i < size_d1 - 1; ++i){
        for( k = 0; k < size_d3 - 1; ++k){
		      fieldB(i,j,k).y() -= dtp * ( rdz * (fieldE(i,j,k+1).x() - fieldE(i,j,k).x() )
                                           - rdx * (fieldE(i+1,j,k).z() - fieldE(i,j,k).z() ) );
        }
    }
    k = size_d3 - 1;
    for(i = 0; i < size_d1 - 1; ++i){
        for( j = 0; j < size_d2 - 1; ++j){
            fieldB(i,j,k).z() -= dtp * ( rdx * (fieldE(i+1,j,k).y() - fieldE(i,j,k).y() )
                                           - rdy * (fieldE(i,j+1,k).x() - fieldE(i,j,k).x() ) );
        }
    }
*/
//}


//void solver_FDTD(Field3d& fieldE, Field3d& fieldB, const Field3d& fieldJ, const World& world){
    // int i, j, k;
    // const double rdx = 1. / Dx;
    // const double rdy = 1. / Dy;
    // const double rdz = 1. / Dz;
    // const int size_d1 = fieldE.size().x();
    // const int size_d2 = fieldE.size().y();
    // const int size_d3 = fieldE.size().z();
    // static Array3D<double> Ez0(2,size_d2,size_d3);
    // static Array3D<double> Ey0(2,size_d2,size_d3);
/*
    for (j = 0; j < size_d2; ++j){
        for (k = 0; k < size_d3; ++k){
    		Ez0(0,j,k) = fieldE(1,j,k).z();
    		Ez0(1,j,k) = fieldE(size_d1-2,j,k).z();	
    		Ey0(0,j,k) = fieldE(1,j,k).y();
    		Ey0(1,j,k) = fieldE(size_d1-2,j,k).y();      
        }
    }

    update_fieldsB(fieldE, fieldB);
	
    i = 0;
    for(j = 1; j < size_d2; ++j){
        for(k = 1; k < size_d3; ++k){
                fieldE(i,j,k).x() += Dt * ( rdy * (fieldB(i,j,k).z() - fieldB(i,j-1,k).z() ) 
                                          - rdz * (fieldB(i,j,k).y() - fieldB(i,j,k-1).y() )
                                          - fieldJ(i,j,k).x() );
        }
    }
	
    for(i = 1; i < size_d1; ++i){
        for(j = 1; j < size_d2; ++j){
            for(k = 1; k < size_d3; ++k){
			    fieldE(i,j,k).x() += Dt * ( rdy * (fieldB(i,j,k).z() - fieldB(i,j-1,k).z() ) 
                                          - rdz * (fieldB(i,j,k).y() - fieldB(i,j,k-1).y() )
                                          - fieldJ(i,j,k).x() );
			    fieldE(i,j,k).y() += Dt * ( rdz * (fieldB(i,j,k).x() - fieldB(i,j,k-1).x() ) 
                                          - rdx * (fieldB(i,j,k).z() - fieldB(i-1,j,k).z() ) 
                                          - fieldJ(i,j,k).y() );
			    fieldE(i,j,k).z() += Dt * ( rdx * (fieldB(i,j,k).y() - fieldB(i-1,j,k).y() ) 
                                          - rdy * (fieldB(i,j,k).x() - fieldB(i,j-1,k).x() ) 
                                          - fieldJ(i,j,k).z() );
		  }
        }
    }
    exchange_fieldsE(fieldE,world.MPIconf);
	
    if(world.region.boundType[0].x() == OPEN ){
        for(j = 0; j < size_d2; ++j){
            for(k = 0; k < size_d3; ++k){
		      double Kabc = (Dt - Dx) / (Dt + Dx);
		      fieldE(0,j,k).z() = Ez0(0,j,k) + Kabc * (fieldE(1,j,k).z() - fieldE(0,j,k).z() );
		      fieldE(0,j,k).y() = Ey0(0,j,k) + Kabc * (fieldE(1,j,k).y() - fieldE(0,j,k).y() );
	  	    }
        }
    }
    
    if(world.region.boundType[1].x() == OPEN){
        for(j = 0; j < size_d2; ++j){
            for(k = 0; k < size_d3; ++k){
	           double Kabc = (Dt - Dx) / (Dt + Dx);
	           fieldE(size_d1-1,j,k).z() = Ez0(1,j,k) + Kabc * (fieldE(size_d1-2,j,k).z() - fieldE(size_d1 - 1,j,k).z() );
	           fieldE(size_d1-1,j,k).y() = Ey0(1,j,k) + Kabc * (fieldE(size_d1-2,j,k).y() - fieldE(size_d1 - 1,j,k).y() );
	       }
        }

    }

   update_fieldsB(fieldE, fieldB);
			*/
//} 

