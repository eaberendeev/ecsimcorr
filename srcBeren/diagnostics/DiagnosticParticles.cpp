#include "Diagnostic.h"

void Writer::write_particles2D(int timestep){ 
  char filename[100];

    for( auto& sp : _species){
    	auto size_x = sp.densityOnGrid.size().x();// - ADD_NODES;
    	auto size_y = sp.densityOnGrid.size().y();
    	auto size_z = sp.densityOnGrid.size().z();
   		Array2D<double> densityPlaneZ(size_x,size_y);
   		Array2D<double> densityPlaneY(size_x,size_z);
   		Array2D<double> densityPlaneX(size_y,size_z);
	    
		int k = size_z/2;
	    for( auto i = 0; i < size_x; i++ ){
	      for( auto j = 0; j < size_y; j++ ){
	            densityPlaneZ(i,j) = sp.densityOnGrid(i,j,k);
	      }
	    }
	   	int j = size_y/2;
	    for( auto i = 0; i < size_x; i++ ){
	        for( auto k = 0; k < size_z; k++ ){
	            densityPlaneY(i,k) = sp.densityOnGrid(i,j,k);
	        }
	    }
		int i = size_x/2;
	    for( auto j = 0; j < size_y; j++ ){
	      for( auto k = 0; k < size_z; k++ ){
	            densityPlaneX(j,k) = sp.densityOnGrid(i,j,k);
	      }
	    }

      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//DensPlaneZ"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneZ, size_x, size_y, filename);
      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//DensPlaneY"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneY, size_x, size_z, filename);
      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//DensPlaneX"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneX, size_y, size_z, filename);


	  densityPlaneZ.clear();
	  size_z -= 3;
	  for( auto i = 0; i < size_x; i++ ){
	      for( auto j = 0; j < size_y; j++ ){
			for( auto k = 0; k < size_z; k++ ){
	            densityPlaneZ(i,j) += float(sp.densityOnGrid(i,j,k)) / size_z;
			}
		  }
	  }

      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//DensPlaneAvgZ"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneZ, size_x, size_y, filename);

      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//CurrentPlaneAvgZ"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
	  write_field2D_AvgPlaneZ(sp.currentOnGrid, filename);
		
	sp.get_Pr();
	
		densityPlaneZ.clear();
		for( auto i = 0; i < size_x; i++ ){
			for( auto j = 0; j < size_y; j++ ){
					densityPlaneZ(i,j) += float(sp.Pxx(i,j)) ;
			}
		}
      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//PxxPlaneAvgZ"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneZ, size_x, size_y, filename);

		densityPlaneZ.clear();
		for( auto i = 0; i < size_x; i++ ){
			for( auto j = 0; j < size_y; j++ ){
					densityPlaneZ(i,j) += float(sp.Pyy(i,j)) ;
			}
		}
      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//PyyPlaneAvgZ"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneZ, size_x, size_y, filename);

	    // for( auto i = 0; i < size_x; i++ ){
	    //   for( auto j = 0; j < size_y; j++ ){
	    //   	for( auto k = 0; k < size_z; k++ ){
		// 		if(fabs(sp.kineticEPredict(i,j,k) ) > 1.e-22){
	    //         	sp.densityOnGrid(i,j,k) = sqrt(1 + sp.Je(i,j,k)/sp.kineticEPredict(i,j,k));
		// 		} else {
		// 			sp.densityOnGrid(i,j,k) = 0.0;
		// 		}

	    //   	}
		//   }
	    // }
	    
	// 	int k1 = size_z/2;
	//     for( auto i = 0; i < size_x; i++ ){
	//       for( auto j = 0; j < size_y; j++ ){
	//             densityPlaneZ(i,j) = sp.densityOnGrid(i,j,k1);
	//       }
	//     }
	//    	int j1 = size_y/2;
	//     for( auto i = 0; i < size_x; i++ ){
	//         for( auto k = 0; k < size_z; k++ ){
	//             densityPlaneY(i,k) = sp.densityOnGrid(i,j1,k);
	//         }
	//     }
	// 	int i1 = size_x/2;
	//     for( auto j = 0; j < size_y; j++ ){
	//       for( auto k = 0; k < size_z; k++ ){
	//             densityPlaneX(j,k) = sp.densityOnGrid(i1,j,k);
	//       }
	//     }

    //   sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//LambdaPlaneZ"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
    //   write_array2D(densityPlaneZ, size_x, size_y, filename);
    //   sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//LambdaPlaneY"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
    //   write_array2D(densityPlaneY, size_x, size_z, filename);
    //   sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//LambdaPlaneX"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
    //   write_array2D(densityPlaneX, size_y, size_z, filename);



    //   sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//Phase2D"+"%04d").c_str(),timestep / TimeStepDelayDiag2D);
    //   write_array2D(sp.phaseOnGrid, sp.phaseOnGrid.size().x(), sp.phaseOnGrid.size().y(),filename);

    }
  
}
