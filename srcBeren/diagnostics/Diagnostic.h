#ifndef DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#include "World.h"
#include "Mesh.h"
#include "Read.h"
#include "ParticlesArray.h"
#include "Timer.h"
#include "util.h"
void write_array2D(const Array2D<double>& data, int size1, int size2, const char* filename);
void write_field2D_AvgPlaneZ(const Field3d& field, const char* filenameCh);

struct DiagDataParams{
    std::vector<double> radiationDiagRadiuses;
    std::vector<double> zondCoordsCircleX;
    std::vector<double> sliceFieldsPlaneX, sliceFieldsPlaneY, sliceFieldsPlaneZ;
    std::vector<double> sliceRadiationPlaneX, sliceRadiationPlaneY, sliceRadiationPlaneZ;
    std::vector<double3> zondCoords;
    std::vector<double3> zondCoordsLineX,zondCoordsLineY,zondCoordsLineZ;
    std::vector<int> outTime3D;

};


struct DiagData{

    
    DiagData(const Region &region){
        //powerRadLine(region ),
          //                  powerRad1D(region.numNodes.x() ){

        std::vector< std::vector<std::string> > vecStringParams;
        read_params_to_string("Diagnostics","./Diagnostics.cfg",vecStringParams);
        for (const auto& line: vecStringParams[0]){
            set_params_from_string(line);
        }
        for( const auto& elem : params.sliceRadiationPlaneX){
            std::cout  << "sliceRadiationPlaneX created for coord " << elem <<"\n";
            powerRadPlaneX.emplace_back( Array2D<double>(region.numNodes.y(), region.numNodes.z() ) );
        }
        for( const auto& elem : params.sliceRadiationPlaneY){
            std::cout  << "sliceRadiationPlaneY created for coord " << elem <<"\n";
            powerRadPlaneY.emplace_back( Array2D<double>(region.numNodes.x(), region.numNodes.z() ) );
        }
        for( const auto& elem : params.sliceRadiationPlaneZ){
            std::cout  << "sliceRadiationPlaneZ created for coord " << elem <<"\n";
            powerRadPlaneZ.emplace_back( Array2D<double>(region.numNodes.x(), region.numNodes.y() ) );
        }        
       Reset() ;

    };
    void set_params_from_string(const std::string& line);

    void calc_energy(Mesh& mesh, const std::vector<ParticlesArray>& species,
                     const ParametersMap& parameters);

    void Reset(){
        
        for(auto& data : powerRadPlaneX){
            data.set_zero();
        }
        for(auto& data : powerRadPlaneY){
            data.set_zero();
        }
        for(auto& data : powerRadPlaneZ){
            data.set_zero();
        }        
    };

    std::vector< Array2D<double> > powerRadPlaneX;
    std::vector< Array2D<double> > powerRadPlaneY;
    std::vector< Array2D<double> > powerRadPlaneZ;

    std::map<std::string,double> energyParticlesKinetic;
    std::map<std::string, double> energyParticlesInject;
    std::map<std::string, double> energyParticlesLost;
    std::map<std::string, double> energy;

    double energyFieldE, energyFieldB, energyFieldBFull;
    double diffEB;
    DiagDataParams params;

};

struct Writer{
protected:
    const World &_world;
    Mesh &_mesh;
    std::vector<ParticlesArray> &_species;
    Domain& _domain;
    ParametersMap& _parameters;

   public:
    FILE *fDiagEnergies;
    DiagData diagData;

    Writer(const World& world, Mesh& mesh,
           std::vector<ParticlesArray>& species, Domain &domain, ParametersMap &parameters);

    void output(double diffV, ParametersMap& parameters, int timestep);
    ~Writer(){
            fclose(fDiagEnergies);  
    } 

    void write_particles2D(int timestep);
    void write_particles3D(int timestep);
    void write_energies(double diffV, const ParametersMap& parameters,
                        int timestep);
    void write_radiation_line(int timestep);
    void write_radiation_avg_line(int timestep);
    
    void write_radiation_circle2D(int timestep);
    void write_radiation_planes(int timestep);
    void write_array2D_planeX(const Array3D<double>& data, double coordX, const std::string& fname, const int& timestep);
    void write_array2D_planeY(const Array3D<double>& data, double coordY, const std::string& fname, const int& timestep);
    void write_array2D_planeZ(const Array3D<double>& data, double coordZ, const std::string& fname, const int& timestep);

    void diag_zond(int timestep);
    void write_fields_lineX(const Field3d& fieldE, const Field3d& fieldB, int timestep);
    void write_fields_lineY(const Field3d& fieldE, const Field3d& fieldB, int timestep);
    void write_fields_lineZ(const Field3d& fieldE, const Field3d& fieldB, int timestep);
    void write_fields3D(const Field3d& fieldE, const Field3d& fieldB,  const int& timestep);
    void write_fields2D_planeX(const Field3d& fieldE, const Field3d& fieldB, double coord, const int& timestep);
    void write_fields2D_planeY(const Field3d& fieldE, const Field3d& fieldB, double coord, const int& timestep);
    void write_fields2D_planeZ(const Field3d& fieldE, const Field3d& fieldB, double coord, const int& timestep);
    void write_fields2D_AvgPlaneZ(const Field3d& fieldE, const Field3d& fieldB, const int& timestep);
    void write_fieldsJ2D_AvgPlaneZ(const Field3d& fieldE, const int& timestep);
    void write_fields2D_planeX(const Field3d& field, double coordX, const std::string& fname, const int& timestep);
    void write_fields2D_planeY(const Field3d& field, double coordY, const std::string& fname, const int& timestep);
    void write_fields2D_planeZ(const Field3d& field, double coordZ, const std::string& fname, const int& timestep);
    void write_fields_circle( int timestep);
    void write_fields2D_circle(const Array2D<double3>& fieldE, const Array2D<double3>& fieldB, int series, const int& timestep);

    void write_fields2D(const Array2D<double3>& fieldE, const Array2D<double3>& fieldB, const int& timestep);
    void write_array2D_planeZ_avg(const Array3D<double>& data,
                                          const std::string& fname,
                                          const int& timestep);
    void write_fields2D_AvgPlaneZ(const Field3d& field,
                                          const std::string& fname,
                                          const int& timestep);
};
#endif 	
