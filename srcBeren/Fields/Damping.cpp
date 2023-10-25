#include "World.h"
#include "Damping.h"
#include "bounds.h"

void Damping_Func(double& source, int i, int maxi, double& energyDamp){
    double koeff=0.8;
    double a, damp;

    a = (1.0-koeff)/(maxi*maxi);
    damp = a * i * i + koeff;

    energyDamp += 0.5*source*source*(1.0 - damp*damp);
    source *= damp;	
} 

double damping_fields(Field3d& fieldE, Field3d& fieldB,const Region& domain){
    int i,j,k,i1,j1,k1;
    double energyDamp, DampCells;
    int max_indx = domain.numNodes.x()-1;
    int max_indy = domain.numNodes.y()-1;
    int max_indz = domain.numNodes.z()-1;

	energyDamp = 0.;

    if (true || !isPeriodicX) { 
      	for(i = 0; i < domain.dampCells[0].x(); i++){
	  		for (j = 0; j <= max_indy ; j++){
		  		for (k = 0; k <= max_indz ; k++){
					
					i1 = i;
					DampCells = domain.dampCells[0].x();
					
					Damping_Func(fieldE(i,j,k,0), i1, DampCells, energyDamp);					
				
					Damping_Func(fieldE(i,j,k,1), i1, DampCells, energyDamp);
								
					Damping_Func(fieldE(i,j,k,2), i1, DampCells, energyDamp);
								
					Damping_Func(fieldB(i,j,k,2), i1, DampCells, energyDamp);

					Damping_Func(fieldB(i,j,k,1), i1, DampCells, energyDamp);				
					
					Damping_Func(fieldB(i,j,k,0), i1, DampCells, energyDamp);
		    	}
	    	}
		}

      	for(i = max_indx; i > max_indx - domain.dampCells[1].x(); i--){
	  	    for (j = 0; j <= max_indy ; j++){
		  	    for (k = 0; k <= max_indz ; k++){

					i1 = -(i - max_indx);
					DampCells = domain.dampCells[1].x();
					
					Damping_Func(fieldE(i,j,k,0), i1, DampCells, energyDamp);							

					Damping_Func(fieldE(i,j,k,1), i1, DampCells, energyDamp);
								
					Damping_Func(fieldE(i,j,k,2), i1, DampCells, energyDamp);
								
					Damping_Func(fieldB(i,j,k,2), i1, DampCells, energyDamp);

					Damping_Func(fieldB(i,j,k,1), i1, DampCells, energyDamp);				
					
					Damping_Func(fieldB(i,j,k,0), i1, DampCells, energyDamp);
		    	}
	    	}
		}
    }

	if (true || !isPeriodicY) { 

	    for( i = 0; i <= max_indx; ++i){			
			for( j = max_indy; j > max_indy - domain.dampCells[1].y(); --j){
				for( k = 0; k <= max_indz; ++k){
								
				    j1 = - j + max_indy;
				    DampCells = domain.dampCells[1].y();

				    Damping_Func(fieldE(i,j,k,0), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,1), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,2), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k,2), j1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,1), j1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,0), j1, DampCells, energyDamp);

				}
			}
			for( j = 0; j< domain.dampCells[0].y(); ++j){
				for( k = 0; k <= max_indz; ++k){
								
				    j1 = j;//- j + domain.numCells_d2+ADD_NODES;
				    DampCells = domain.dampCells[0].y();

				    Damping_Func(fieldE(i,j,k,0), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,1), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,2), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k,2), j1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,1), j1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,0), j1, DampCells, energyDamp);
				}
			}		
	    }
  	}

	if (true || !isPeriodicZ) { 

	    for( i = 0; i <= max_indx; ++i){			
			for( j = 0; j <= max_indy; ++j){
				for( k = max_indz; k > max_indz - domain.dampCells[1].z(); --k){
								
				    k1 = - k + max_indz;
				    DampCells = domain.dampCells[1].z();

				    Damping_Func(fieldE(i,j,k,0), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,1), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,2), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k,2), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,1), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,0), k1, DampCells, energyDamp);
				}
			
			for( k = 0; k< domain.dampCells[0].z(); ++k){
								
				    k1 = k;//- j + domain.numCells_d2+ADD_NODES;
				    DampCells = domain.dampCells[0].z();

				    Damping_Func(fieldE(i,j,k,0), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,1), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k,2), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k,2), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,1), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k,0), k1, DampCells, energyDamp);
				}
			}		
	    }
  	}
	return energyDamp;
}
