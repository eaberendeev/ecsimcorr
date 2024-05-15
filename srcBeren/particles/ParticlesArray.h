#ifndef PARTICLES_ARRAY_H_
#define PARTICLES_ARRAY_H_
#include <assert.h>

#include <functional>

#include "Mesh.h"
#include "Particle.h"
#include "World.h"
#include "containers.h"
#include "random_generator.h"

typedef Eigen::Triplet<double> Trip;


class ParticlesArray{

public:
    Array3D< std::vector<Particle> > particlesData;
    Array3D<int> countInCell;

    Array3D<double> densityOnGrid;
    Array2D<double> phaseOnGrid;
    Field3d currentOnGrid;
    Array3D<double> Pxx;
    Array3D<double> Pyy;
    Array3D<double> Pzz;

    const double charge;
    const double density;
    double phasePXmin, phasePXmax;
    double kineticEnergy;
    double injectionEnergy;
    double lostEnergy;
    const std::string _name;
    const double temperature;
    const int NumPartPerCell;
    void delete_bounds();
    void add_particle(const Particle &particle);
    void save_init_coord();
    void save_init_velocity();
    void save_init_coord_and_velocity();
    void calc_Esirkepov_current(const double dt, Field3d& fieldJ) const;
    void update_cells(const Domain& domain);
    void set_particles();

    const std::string& name() const noexcept { return _name; }

    void delete_particle_runtime(int ix, int iy, int iz, int ip){
        countInCell(ix,iy,iz)--;
        int old_count = countInCell(ix,iy,iz);
        particlesData(ix,iy,iz)[ip] = particlesData(ix,iy,iz)[old_count];
        int lastParticle = particlesData(ix,iy,iz).size()-1;
        if(old_count == lastParticle ){
             particlesData(ix,iy,iz).pop_back();
        }
        else{
            particlesData(ix,iy,iz)[old_count] = particlesData(ix,iy,iz)[lastParticle];
            particlesData(ix,iy,iz).pop_back();
        }
    }
    void update_count_in_cell(){
        for (auto i = 0;i< countInCell.size().x();i++){
            for (auto j = 0; j< countInCell.size().y();j++){
                for (auto k = 0; k< countInCell.size().z();k++){
                    countInCell(i,j,k) = particlesData(i,j,k).size();
                }
            }
        }
    }

    void add_uniform_line(int numParts, double3 x0, double3 size);

    void set_test_particles();
    void set_distribution_density(std::function<double(double3 )> set_density);
    void set_smooth_mass();
    ParticlesArray(const ParametersMap& particlesParams,
                   const ParametersMap& parameters, const Domain& domain);
    void set_params_from_string(const std::string& line);
    void density_on_grid_update();
    void phase_on_grid_update(const Domain& domain);
    void inject(int timestep);
    void update(Mesh& mesh,int timestep);
    double mass() const{
          return _mass;
    }

    double mpw() const{
          return _mpw;
    }    
    std::vector<Particle>& operator() (int i) {
        return particlesData(i);
    }

    const std::vector<Particle>& operator() (int i) const{
        return particlesData(i);
    }
    int size() const{
        return particlesData.size().x()*particlesData.size().y()*particlesData.size().z();
    }
    double get_kinetic_energy() const;
    double3 get_kinetic_energy_component() const;
    double get_kinetic_energy(int dim) const;
    double get_kinetic_energy(int dim1, int dim2) const;
    double get_inject_energy() const { return injectionEnergy; }
    void glue_density_bound();
    void move(double dt);
    void prepare(int timestep);
    void correctv(Mesh& mesh, const Domain& domain, const double dt);
    void correctv_component(Mesh& mesh, const Domain& domain, const double dt);
    void predict_velocity(const Mesh& mesh, const Domain& domain,
                          const double dt);
    void move_and_calc_current(const double dt, Field3d& fieldJ);
    void move_and_calc_current(const double dt);
    void predict_current(const Field3d& fieldB, Field3d& fieldJ,
                         const Domain& domain, const double dt);
    void get_L(Mesh& mesh, const Domain& domain, const double dt);
    bool particle_boundaries(double3& coord, const Domain& domain);
    void get_P();
    void get_Pr();
    int get_total_num_of_particles() {
        int count = 0;
        for (auto k = 0; k < size(); ++k) {
            count += particlesData(k).size();
        }

        return count;
    }
    double add_uniform_cilinderZ(const int numParts, const double3& temperature,
                                 const double3& c, const double r0,
                                 const double z0,
                                 ThreadRandomGenerator& randGenSpace,
                                 ThreadRandomGenerator& randGenPulse);
    double add_uniform_rectangle(const int numParts, const double3& temperature,
                                 const double3& startsCoord,
                                 const double3& rectSize,
                                 ThreadRandomGenerator& randGenSpace,
                                 ThreadRandomGenerator& randGenPulse);

   protected:
    double _mass;
    double _mpw; /*macroparticle weight*/
    double xCellSize;
    double yCellSize;
    double zCellSize;
    int xCellCount;
    int yCellCount;
    int zCellCount;
};



double PulseFromKev(double kev, double mass);

int get_num_of_type_particles(const std::vector<ParticlesArray> &Particles, const std::string& ParticlesType);
/// Ionization particles = electron(particles_e) +  ion (particles_i) 
void collision(const Mesh &mesh, std::vector<ParticlesArray> &Particles, int timestep);
void reserve_Lmat(Mesh& mesh, ParticlesArray &sp);

#endif 
