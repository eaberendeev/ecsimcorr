#include "Diagnostic.h"

// Устанавливаем значения основных параметров в переменной Params в соответствии с содержимым сроки
void DiagData::set_params_from_string(const std::string& line){
    std::vector<std::string> strvec;

    strvec = split(line, ' ');

    if(strvec[0]=="radiationDiagRadiuses") {
      for(uint i = 1; i < strvec.size();i++ ){
        this->params.radiationDiagRadiuses.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceRadiationPlaneX") {
      for(uint i = 1; i < strvec.size(); i++ ){
        this->params.sliceRadiationPlaneX.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceRadiationPlaneY") {
      for(uint i = 1; i < strvec.size(); i++ ){
        this->params.sliceRadiationPlaneY.push_back( stod(strvec[i] ) );
      }
    }

    if(strvec[0]=="sliceRadiationPlaneZ") {
      for(uint i = 1; i < strvec.size() ; i++){
        this->params.sliceRadiationPlaneZ.push_back( stod(strvec[i] ) );
      }
    }

    if(strvec[0]=="sliceFieldsPlaneX") {
      for(uint i = 1; i < strvec.size() ; i++){
        this->params.sliceFieldsPlaneX.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceFieldsPlaneY") {
      for(uint i = 1; i < strvec.size(); i++ ){
        this->params.sliceFieldsPlaneY.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceFieldsPlaneZ") {
      for(uint i = 1; i < strvec.size(); i++){
        this->params.sliceFieldsPlaneZ.push_back( stod(strvec[i] ) );
      }
    }   
    if(strvec[0]=="zondCoords") {
      for(uint i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoords.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }
    if(strvec[0]=="zondCoordsLineX") {
      for(uint i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoordsLineX.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }          
    if(strvec[0]=="zondCoordsLineY") {
      for(uint i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoordsLineY.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }   
    if(strvec[0]=="zondCoordsLineZ") {
      for(uint i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoordsLineZ.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }       

    if(strvec[0]=="outTime3D") {
      for(uint i = 1; i < strvec.size(); i++ ){
        this->params.outTime3D.push_back( stol(strvec[i] ) );
      }
    }
    

}


void Writer::output(double diffV,  int timestep){

////////////////////////////////////////////////
////////////////////////////////////////////////
//////////////// WRITE DATA ////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////    


/////// 2D DATA //////////////////////////////
    if (timestep % TimeStepDelayDiag2D == 0){

        write_particles2D(timestep);
        for (auto coordX : diagData.params.sliceFieldsPlaneX){
            write_fields2D_planeX(_mesh.fieldEn, _mesh.fieldB, coordX, timestep);
        }
        for (auto coordY : diagData.params.sliceFieldsPlaneY){
            write_fields2D_planeY(_mesh.fieldEn, _mesh.fieldB, coordY, timestep);
        }     
        for (auto coordZ : diagData.params.sliceFieldsPlaneZ){
            write_fields2D_planeZ(_mesh.fieldEn, _mesh.fieldB, coordZ, timestep);
        }
        //write_fields2D_AvgPlaneZ(_mesh.fieldEn, _mesh.fieldB, timestep);

    }

        diagData.calc_energy(_mesh,_species);

        write_energies(diffV, timestep);

}


void make_folders(){
    mkdir(".//Fields", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Fields//Diag3D", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Fields//Diag2D", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Fields//Diag1D", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Recovery", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Recovery//Fields", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Recovery//Particles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Anime", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Particles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Performance", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::cout << "Create folders for output...\n";
}

Writer::Writer(const World &world, const Mesh &mesh,std::vector<ParticlesArray> &species) : 
    _world(world),_mesh(mesh),_species(species),diagData(world.region) {
  
  make_folders(); 
  for( const auto &sp : _species){
    mkdir((".//Particles//" + sp.name).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((".//Particles//" + sp.name+"//Diag1D").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((".//Particles//" + sp.name+"//Diag2D").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((".//Particles//" + sp.name+"//Diag3D").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((".//Recovery//Particles//" + sp.name).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  fDiagEnergies = fopen("Energies.dat", "w");
  
}

void write_array2D(const Array2D<double>& data, int size1, int size2, const char* filenameCh){
    
    std::string filename(filenameCh);
    std::ofstream fData2D(filename, std::ios::out | std::ios::binary);

    float info;
    
    Array2D<float> floatData(size1,size2);
    
    for( auto i = 0; i < size1; ++i ){
       for( auto j = 0; j < size2; ++j ){
          floatData(i,j) = float(data(i,j));
      }
    }

          info = float(size1);
          fData2D.write((char*) &info, sizeof(info));
          info = float(size2);
          fData2D.write((char*) &info, sizeof(info));
        
     fData2D.write((char*) &floatData.data(0), size1*size2 * sizeof(floatData.data(0)) );
 
}


void write_field2D_AvgPlaneZ(const Field3d& field, const char* filenameCh){
    float info;    
    int indx;

    int size_x = field.size().x(); // - ADD_NODES;
    int size_y = field.size().y(); // - ADD_NODES;
    int size_z = field.size().z();
    int size1 = size_x;
    int size2 = size_y;


    float* floatData[3]; 

    for(auto i = 0; i<3; i++){
        floatData[i] = new float[size1*size2];
    }
      
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
            for( auto k = 0; k<3; k++){
            indx = i*size_y + j;
            floatData[k][indx] = 0;
          }
        }
      }

    std::string filename(filenameCh);
    std::ofstream fField2D(filename, std::ios::out | std::ios::binary);    
    
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
        for( auto k = 0; k < size_z; k++ ){
              indx = i*size_y + j;
              floatData[0][indx] += float(field(i,j,k,0) ) / size_z;
              floatData[1][indx] += float(field(i,j,k,1) ) / size_z;
              floatData[2][indx] += float(field(i,j,k,2) ) / size_z;
        }
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
