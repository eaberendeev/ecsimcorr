#ifndef PARTICLES_ARRAY_H_
#define PARTICLES_ARRAY_H_
#include "Particle.h"
#include "World.h"
#include "Vec.h"
#include "Mesh.h"
#include <functional>
#include <assert.h>
#include "random_generator.h"

typedef Eigen::Triplet<double> Trip;

struct ParticlesOption{
    int boundResumption;
    int sourceType;
    int smoothMass;
    double smoothMassSize;
    double smoothMassMax;
    double sourceAngle; 
};


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

    double charge;
    double density;
    double phasePXmin, phasePXmax;
    double pot_I;
    double pot_k;
    double kineticEnergy;
    double injectionEnergy;
    double lostEnergy;
    std::string name;
    double temperature;
    double velocity;
    double widthY,widthZ;
    double focus;
    static int counter;
    void delete_bounds();
    ParticlesOption option;
    std::string initDist;
    std::vector<double> distParams;
    void add_particle(const Particle &particle);
    void save_init_coord();
    void calc_Esirkepov_current(const double dt, Field3d& fieldJ) const;
    void update_cells();

    void write_particles_to_recovery(int timestep);
    void read_particles_from_recovery();

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
    ParticlesArray(const std::vector<std::string>& vecStringParams, World& world);
    void set_params_from_string(const std::string& line);
    void density_on_grid_update();
    void phase_on_grid_update();
    void set_distribution();
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
    void correctv(Mesh& mesh);
    void correctv_component(Mesh& mesh);
    void predict_velocity(const Mesh& mesh);
    void move_and_calc_current(const double dt, Field3d& fieldJ);
    void move_and_calc_current(const double dt);
    void predict_current(const Field3d& fieldB, Field3d& fieldJ);
    void get_L( Mesh& mesh);
    bool particle_boundaries(double3& coord);
    void get_P();
    void get_Pr();

    void set_space_distribution();
    void set_pulse_distribution(ThreadRandomGenerator& randGen);
    void set_uniform_circle(int3 start, int3 end);
    void set_strict_uniform(int3 start, int3 end);
    void set_uniform(int3 start, int3 end, RandomGenerator& randGen);
    void inject_particles(const int timestep);
    double add_uniform_cilinder(int numParts, double r0, double z0, double3 c,
                                ThreadRandomGenerator& randGenSpace,
                                ThreadRandomGenerator& randGenPulse);
    double add_uniform_line(int numParts, double3 r, double3 sizeL,
                            ThreadRandomGenerator& randGenSpace,
                            ThreadRandomGenerator& randGenPulse);

   protected:
    World &_world;
    double _mass;
    double _mpw; /*macroparticle weight*/
};



double PulseFromKev(double kev, double mass);

int get_num_of_type_particles(const std::vector<ParticlesArray> &Particles, const std::string& ParticlesType);
/// Ionization particles = electron(particles_e) +  ion (particles_i) 
void collision(const Mesh &mesh, const World& world ,std::vector<ParticlesArray> &Particles,int timestep);
void reserve_Lmat(Mesh& mesh, ParticlesArray &sp);

#endif 
