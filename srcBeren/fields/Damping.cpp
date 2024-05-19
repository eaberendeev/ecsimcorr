#include "World.h"
#include "Damping.h"
#include "bounds.h"

/*void SetEnergyDampLineNull(const Region& domain){
    for( auto i = 0; i < NumCellsR_glob + 2 * SHAPE_SIZE - 1; ++i){
	energy.eDampLineB_Z1[i] = 0.;
	energy.eDampLineE_Z1[i] = 0.;
	energy.eDampLineB_Z2[i] = 0.;
	energy.eDampLineE_Z2[i] = 0.;
    }
    for( auto i = 0; i < NumCellsZ_glob + 2 * SHAPE_SIZE - 1; ++i){
	energy.eDampLineB_R[i] = 0.;
	energy.eDampLineE_R[i] = 0.;
    }
}*/

void Damping_Func(double& source, int i, int maxi, double& energyDamp){
    double koeff=0.8;
    double a, damp;

    a = (1.0-koeff)/(maxi*maxi);
    damp = a * i * i + koeff;

    energyDamp += 0.5*source*source*(1.0 - damp*damp);
    source *= damp;	
} 
void Damping_Func(double& source, double i, double maxi, double& energyDamp){
    double koeff=0.8;
    double a, damp;
	if(i >= maxi) return;

    a = (1.0-koeff)/(maxi*maxi);
    damp = a * i * i + koeff;

    energyDamp += 0.5*source*source*(1.0 - damp*damp);
    source *= damp;	
} 
double damping_fields(Field3d& fieldE, Field3d& fieldB,const Region& domain){
    int i,j,k; //,i1,j1,k1;
    double energyDamp, DampCells;
    //int max_indx = domain.numNodes.x()-1;
    //int max_indy = domain.numNodes.y()-1;
    //int max_indz = domain.numNodes.z()-1;
	DampCells=10*Dx;
	double center = 0.5*Dx*(fieldE.size().x()-1);
	int centeri = (fieldE.size().x()-1)/2;

	energyDamp = 0.;

        for(i = 0; i < fieldE.size().x(); i++){
                        for (j = 0; j < fieldE.size().y() ; j++){
                                for (k = 0; k < fieldE.size().z() ; k++){
	double r = sqrt( Dx*Dx*( (i-centeri)*(i-centeri)+(j-centeri)*(j-centeri) ) );
if(r <= center - DampCells) {
continue;
} else if( r>= center){
fieldE(i,j,k,0) *=0.8;
fieldE(i,j,k,1) *=0.8;
fieldE(i,j,k,2) *=0.8;
fieldB(i,j,k,0) *=0.8;
fieldB(i,j,k,1) *=0.8;
fieldB(i,j,k,2) *=0.8;
} else{
Damping_Func(fieldE(i,j,k,0), center - r, DampCells, energyDamp);
Damping_Func(fieldE(i,j,k,1), center - r, DampCells, energyDamp);
Damping_Func(fieldE(i,j,k,2), center-r, DampCells, energyDamp);
Damping_Func(fieldB(i,j,k,0), center-r, DampCells, energyDamp);
Damping_Func(fieldB(i,j,k,1), center-r, DampCells, energyDamp);
Damping_Func(fieldB(i,j,k,2), center-r, DampCells, energyDamp);
}					

}
}
}
/*

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

	if (!isPeriodicZ) { 

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
  	}*/
	return energyDamp;
}
