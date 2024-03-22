#include "Diagnostic.h"
#include <vector>

void Writer::diag_zond(int timestep){
  
  char filename[100];
  static FILE *fZond; 
  uint n;
  double xx,yy,zz;  
  double3 E,B,r;


  if (timestep == StartTimeStep){
    sprintf(filename, "./Fields/Zond%04d.dat",0 );
    fZond = fopen(filename, "w");
    fprintf(fZond, "%s ", "## timestep ");
    for (n = 0; n < diagData.params.zondCoords.size(); ++n){
      
      if( ! _world.region.in_region(diagData.params.zondCoords[n].x() ) ) continue;
      xx = diagData.params.zondCoords[n].x();
      yy = diagData.params.zondCoords[n].y();
      zz = diagData.params.zondCoords[n].z();
      fprintf(fZond, "%s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s ", 
          "Ex(",xx,",",yy,",",zz,")","Ey(",xx,",",yy,",",zz,")","Ez(",xx,",",yy,",",zz,")",
          "Bx(",xx,",",yy,",",zz,")","By(",xx,",",yy,",",zz,")","Bz(",xx,",",yy,",",zz,")"  );
    }
    fprintf(fZond, "\n");
  }
  
    fprintf(fZond, "%g ",Dt*timestep);
    
    for (n = 0; n < diagData.params.zondCoords.size(); ++n){
      if( ! _world.region.in_region( diagData.params.zondCoords[n].x() ) ) continue;
      r.x() = diagData.params.zondCoords[n].x() -  _world.region.origin;
      r.y() = diagData.params.zondCoords[n].y();
      r.z() = diagData.params.zondCoords[n].z();

      E = get_fieldE_in_pos(_mesh.fieldE,r);
      B = get_fieldB_in_pos(_mesh.fieldB,r);
      fprintf(fZond, "%g %g %g %g %g %g ",  E.x(), E.y(), E.z(), B.x(), B.y(), B.z() );
    }
    
    fprintf(fZond, "\n");
    if(  timestep % TimeStepDelayDiag1D == 0){
      fflush(fZond);
    }

}


void Writer::write_fields2D_planeX(const Field3d& fieldE, const Field3d& fieldB, double coordX, const int& timestep){

    int globIndex = _mesh.get_node_from_coordX(coordX);
    int index = _world.region.get_index_loc(globIndex);
    
    char filenameCh[100];
    float info;    
    int indx;

    int size_y = fieldE.size().y(); //- ADD_NODES;
    int size_z = fieldE.size().z(); //- ADD_NODES;
    int size1 = size_y;
    int size2 = size_z;

    float* floatData[6];

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    sprintf(filenameCh, ".//Fields//Diag2D//FieldPlaneX_%04d_%04d",globIndex,timestep / TimeStepDelayDiag2D); 
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);
        
     int i = index; 
      for( auto j = 0; j < size_y; j++ ){
          for( auto k = 0; k < size_z; k++ ){
              indx = j*size_z + k;
              floatData[0][indx] = float(fieldE(i,j,k,0) );
              floatData[1][indx] = float(fieldE(i,j,k,1) );
              floatData[2][indx] = float(fieldE(i,j,k,2) );
              floatData[3][indx] = float(fieldB(i,j,k,0) );
              floatData[4][indx] = float(fieldB(i,j,k,1) );
              floatData[5][indx] = float(fieldB(i,j,k,2) );
          }
      }
    
        info = float(size1);
        fField2D.write((char*) &info, sizeof(info));
        info = float(size2);
        fField2D.write((char*) &info, sizeof(info));

        for(auto i = 0; i<6; ++i){
          fField2D.write((char*) floatData[i], size1*size2 * sizeof(float) );
        }
        
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}

void Writer::write_fields2D_AvgPlaneZ(const Field3d& fieldE, const Field3d& fieldB, const int& timestep){
    char filenameCh[100];
    float info;    
    int indx;

    int size_x = fieldE.size().x(); // - ADD_NODES;
    int size_y = fieldE.size().y(); // - ADD_NODES;
    int size_z = fieldE.size().z(); // - ADD_NODES;
    //int size_z = fieldE.size().z();
    int size1 = size_x;
    int size2 = size_y;


    float* floatData[6]; 

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }
      
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
            for( auto k = 0; k<6; k++){
            indx = i*size_y + j;
            floatData[k][indx] = 0;
          }
        }
      }

    sprintf(filenameCh, ".//Fields//Diag2D//FieldAvgPlaneZ_%04d",timestep / TimeStepDelayDiag2D);
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);    
    
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
        for( auto k = 1; k < size_z - 2; k++ ){
              indx = i*size_y + j;
              floatData[0][indx] +=
                  float(fieldE(i, j, k, 0) / (size_z - ADD_NODES));
              floatData[1][indx] +=
                  float(fieldE(i, j, k, 1) / (size_z - ADD_NODES));
              floatData[2][indx] +=
                  float(fieldE(i, j, k, 2) / (size_z - ADD_NODES));
              floatData[3][indx] +=
                  float(fieldB(i, j, k, 0) / (size_z - ADD_NODES));
              floatData[4][indx] +=
                  float(fieldB(i, j, k, 1) / (size_z - ADD_NODES));
              floatData[5][indx] += float(fieldB(i,j,k,2)  / (size_z-ADD_NODES));
        }
      }
    }
    

        info = float(size1);
        fField2D.write((char*) &info, sizeof(info));
        info = float(size2);
        fField2D.write((char*) &info, sizeof(info));
    
    for(auto i = 0; i<6; ++i){
          fField2D.write((char*) floatData[i], size1*size2 * sizeof(float) );
    }
        
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}

void Writer::write_fields2D_planeZ(const Field3d& fieldE, const Field3d& fieldB, double coordZ, const int& timestep){
    char filenameCh[100];
    float info;    

    int indx;

    int size_x = fieldE.size().x(); // - ADD_NODES;
    int size_y = fieldE.size().y(); // - ADD_NODES;
    int size1 = size_x;
    int size2 = size_y;


    float* floatData[6]; 

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    sprintf(filenameCh, ".//Fields//Diag2D//FieldPlaneZ_%04d_%04d",_mesh.get_node_from_coordZ(coordZ),timestep / TimeStepDelayDiag2D);    
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);    

    int k = _mesh.get_node_from_coordZ(coordZ); 
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
            indx = i*size_y + j;
              floatData[0][indx] = float(fieldE(i,j,k,0) );
              floatData[1][indx] = float(fieldE(i,j,k,1) );
              floatData[2][indx] = float(fieldE(i,j,k,2) );
              floatData[3][indx] = float(fieldB(i,j,k,0) );
              floatData[4][indx] = float(fieldB(i,j,k,1) );
              floatData[5][indx] = float(fieldB(i,j,k,2) );
        }
      }
    
        info = float(size1);
        fField2D.write((char*) &info, sizeof(info));
        info = float(size2);
        fField2D.write((char*) &info, sizeof(info));
    
    for(auto i = 0; i<6; ++i){
          fField2D.write((char*) floatData[i], size1*size2 * sizeof(float) );
    }
        
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}


void Writer::write_fields2D_planeY(const Field3d& fieldE, const Field3d& fieldB, double coordY, const int& timestep){
    char filenameCh[100];
    float info;    

    int indx;

    int size_x = fieldE.size().x(); // - ADD_NODES;
    int size_z = fieldE.size().z(); //- ADD_NODES;
    int size1 = size_x;
    int size2 = size_z;

    float* floatData[6]; 

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    sprintf(filenameCh, ".//Fields//Diag2D//FieldPlaneY_%04d_%04d",_mesh.get_node_from_coordY(coordY),timestep / TimeStepDelayDiag2D); 
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);    

     int j = _mesh.get_node_from_coordY(coordY); 
      for( auto i = 0; i < size_x; i++ ){
          for( auto k = 0; k < size_z; k++ ){
            indx = i*size_z + k;
              floatData[0][indx] = float(fieldE(i,j,k,0) );
              floatData[1][indx] = float(fieldE(i,j,k,1) );
              floatData[2][indx] = float(fieldE(i,j,k,2) );
              floatData[3][indx] = float(fieldB(i,j,k,0) );
              floatData[4][indx] = float(fieldB(i,j,k,1) );
              floatData[5][indx] = float(fieldB(i,j,k,2) );
          }
      }

        info = float(size1);
        fField2D.write((char*) &info, sizeof(info));
        info = float(size2);
        fField2D.write((char*) &info, sizeof(info));
    
    for(auto i = 0; i<6; ++i){
          fField2D.write((char*) floatData[i], size1*size2 * sizeof(float) );
    }
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}


void Writer::write_fields2D_planeX(const Field3d& field, double coordX, const std::string& fname, const int& timestep){

    int globIndex = _mesh.get_node_from_coordX(coordX);
    int index = _world.region.get_index_loc(globIndex);
    
    char filenameCh[100];
    float info;    
    int indx;

    int size_y = field.size().y(); //- ADD_NODES;
    int size_z = field.size().z(); //- ADD_NODES;
    int size1 = size_y;
    int size2 = size_z;

    float* floatData[3];

    for(auto i = 0; i<3; i++){
        floatData[i] = new float[size1*size2];
    }

    sprintf(filenameCh, (fname+"PlaneX_%04d_%04d").c_str(),globIndex,timestep / TimeStepDelayDiag2D); 
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);
        
     int i = index; 
      for( auto j = 0; j < size_y; j++ ){
          for( auto k = 0; k < size_z; k++ ){
              indx = j*size_z + k;
              floatData[0][indx] = float(field(i,j,k,0) );
              floatData[1][indx] = float(field(i,j,k,1) );
              floatData[2][indx] = float(field(i,j,k,2) );
          }
      }
    
        info = float(size1);
        fField2D.write((char*) &info, sizeof(info));
        info = float(size2);
        fField2D.write((char*) &info, sizeof(info));

        for(auto i = 0; i<3; ++i){
          fField2D.write((char*) floatData[i], size1*size2 * sizeof(float) );
        }
        
    for(auto i = 0; i<3; i++){
        delete[] floatData[i];
    }

}

void Writer::write_fields2D_planeZ(const Field3d& field, double coordZ, const std::string& fname, const int& timestep){
    char filenameCh[100];
    float info;    

    int indx;

    int size_x = field.size().x(); // - ADD_NODES;
    int size_y = field.size().y(); // - ADD_NODES;
    int size1 = size_x;
    int size2 = size_y;


    float* floatData[3]; 

    for(auto i = 0; i<3; i++){
        floatData[i] = new float[size1*size2];
    }
    int globIndex = _mesh.get_node_from_coordZ(coordZ);

    sprintf(filenameCh, (fname+"PlaneZ_%04d_%04d").c_str(),globIndex,timestep / TimeStepDelayDiag2D); 
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);    

    int k = _mesh.get_node_from_coordZ(coordZ); 
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
            indx = i*size_y + j;
              floatData[0][indx] = float(field(i,j,k,0) );
              floatData[1][indx] = float(field(i,j,k,1) );
              floatData[2][indx] = float(field(i,j,k,2) );
        }
      }
    
        info = float(size1);
        fField2D.write((char*) &info, sizeof(info));
        info = float(size2);
        fField2D.write((char*) &info, sizeof(info));
    
    for(auto i = 0; i<3; ++i){
          fField2D.write((char*) floatData[i], size1*size2 * sizeof(float) );
    }
        
    for(auto i = 0; i<3; i++){
        delete[] floatData[i];
    }

}

void Writer::write_fields2D_AvgPlaneZ(const Field3d& field,
                                   const std::string& fname,
                                   const int& timestep) {
    char filenameCh[100];
    float info;

    int indx;

    int size_x = field.size().x();   // - ADD_NODES;
    int size_y = field.size().y();   // - ADD_NODES;
    int size_z = field.size().z();   // - ADD_NODES;
    int size1 = size_x;
    int size2 = size_y;

    float* floatData[3];

    for (auto i = 0; i < 3; i++) {
        floatData[i] = new float[size1 * size2];
    }

    sprintf(filenameCh, (fname + "PlaneAvgZ_%04d").c_str(),
            timestep / TimeStepDelayDiag2D);
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);
    for (auto i = 0; i < size_x; i++) {
        for (auto j = 0; j < size_y; j++) {
            for (auto k = 0; k < size_z; k++) {
                floatData[0][indx] = 0;
                floatData[1][indx] = 0;
                floatData[2][indx] = 0;
            }
        }
    }

    for (auto i = 0; i < size_x; i++) {
        for (auto j = 0; j < size_y; j++) {
            for (auto k = 1; k < size_z - 2; k++) {
                indx = i * size_y + j;
                floatData[0][indx] += float(field(i, j, k, 0)/(size_z-ADD_NODES));
                floatData[1][indx] += float(field(i, j, k, 1)/(size_z-ADD_NODES));
                floatData[2][indx] += float(field(i, j, k, 2)/(size_z-ADD_NODES));
            }
        }
    }

    info = float(size1);
    fField2D.write((char*) &info, sizeof(info));
    info = float(size2);
    fField2D.write((char*) &info, sizeof(info));

    for (auto i = 0; i < 3; ++i) {
        fField2D.write((char*) floatData[i], size1 * size2 * sizeof(float));
    }

    for (auto i = 0; i < 3; i++) {
        delete[] floatData[i];
    }
}

void Writer::write_fields2D_planeY(const Field3d& field, double coordY, const std::string& fname, const int& timestep){
    char filenameCh[100];
    float info;    

    int indx;

    int size_x = field.size().x(); // - ADD_NODES;
    int size_z = field.size().z(); //- ADD_NODES;
    int size1 = size_x;
    int size2 = size_z;

    float* floatData[3]; 

    for(auto i = 0; i<3; i++){
        floatData[i] = new float[size1*size2];
    }
    int globIndex = _mesh.get_node_from_coordY(coordY);

    sprintf(filenameCh, (fname+"PlaneY_%04d_%04d").c_str(),globIndex,timestep / TimeStepDelayDiag2D); 
    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);    

     int j = _mesh.get_node_from_coordY(coordY); 
      for( auto i = 0; i < size_x; i++ ){
          for( auto k = 0; k < size_z; k++ ){
            indx = i*size_z + k;
              floatData[0][indx] = float(field(i,j,k,0) );
              floatData[1][indx] = float(field(i,j,k,1) );
              floatData[2][indx] = float(field(i,j,k,2) );
          }
      }

        info = float(size1);
        fField2D.write((char*) &info, sizeof(info));
        info = float(size2);
        fField2D.write((char*) &info, sizeof(info));
    
    for(auto i = 0; i<3; ++i){
          fField2D.write((char*) floatData[i], size1*size2 * sizeof(float) );
    }
    
    for(auto i = 0; i<3; i++){
        delete[] floatData[i];
    }

}
