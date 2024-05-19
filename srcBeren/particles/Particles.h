#ifndef PARTICLES_H_
#define PARTICLES_H_
#include "World.h"
#include "Vec.h"
#include "Mesh.h"
#include <functional>
#include <assert.h>

typedef Eigen::Triplet<double> Trip;

struct ParticleSimple{
	double3 coord;
    double3 velocity;
	double3 initCoord;
	double3 initVelocity;

	friend std::ostream& operator<<(std::ostream& out, const ParticleSimple &particle);
    
    void set_global(const Region& domain){
        coord.x() += domain.origin;
    }
    void set_local(const Region& domain){
        coord.x() -= domain.origin;      
    }
    void move(double dt){
        coord += velocity * dt;
    }
};

    inline int pos_ind(int index, int n, int _size1, int _size2, int _size3){
        std::vector<int> dim = {_size1, _size2, _size3};
        int capacity = 1;
        for(unsigned int i = n + 1; i < dim.size(); i++){
            capacity *= dim[i];
        }
        return (index / capacity) % dim[n];
    }

struct ParticleMPW : ParticleSimple{
    double mpw;
    friend std::ostream& operator<<(std::ostream& out, const ParticleMPW &particle);

};
struct ParticleMass : ParticleSimple{
	double mass;
  	friend std::ostream& operator<<(std::ostream& out, const ParticleMass &particle);
};

struct ParticleMassMPW : ParticleSimple{
    double mass,mpw;
    friend std::ostream& operator<<(std::ostream& out, const ParticleMassMPW &particle);
};

#ifdef PARTICLE_MASS
    #ifdef PARTICLE_MPW
        typedef ParticleMassMPW Particle;
    #else
        typedef ParticleMass Particle;
    #endif
#else
    #ifdef PARTICLE_MPW
        typedef ParticleMPW Particle;
    #else
        typedef ParticleSimple Particle;
    #endif
#endif


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
    //Array3D<double> Je;
    //Array3D<double> kineticEPredict;
    Field3d currentOnGrid;
    Array2D<double> Pxx;
    Array2D<double> Pyy;

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
    void correctv2(const Field3d& fieldB,const Field3d& fieldE);

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

    void stencil_Lmat( Mesh& mesh);
    void set_test_particles();
    void set_distribution_density(std::function<double(double3 )> set_density);
    void set_smooth_mass();
    ParticlesArray(const std::vector<std::string>& vecStringParams, World& world);
    void set_params_from_string(const std::string& line);
    void density_on_grid_update();
    //void density_on_grid_update_pic();
    void phase_on_grid_update();
    void set_distribution();
    void inject(int timestep);
    void update(Mesh& mesh,int timestep);
    double mass() const{
          return _mass;
    }
    // double mass(int k) const{
    //     #ifdef PARTICLE_MASS
    //       return particlesData(k).mass;
    //     #else
    //       return _mass;
    //     #endif
    // }
    double mpw() const{
        // #ifdef PARTICLE_MPW
        //   return particlesData(k).mpw;
        // #else
          return _mpw;
        //#endif
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
    double get_inject_energy() const{
        return injectionEnergy;
    }
    void glue_density_bound();
    //void move(Mesh& mesh,int timestep);
    void move(double dt);
    void prepare(int timestep);
    //void correctv(const double dt, Field3d& fieldE, Field3d& fieldEn, Field3d& fieldEp, Field3d& fieldJp);
    void correctv(Mesh& mesh);
    void predict_velocity(const Mesh& mesh);
    void move_and_calc_current(const double dt, Field3d& fieldJ);
    void move_and_calc_current(const double dt);
    void predict_current(const Field3d& fieldB, Field3d& fieldJ);
    void predict_current2( const Field3d& fieldB, Field3d& fieldJ);
    void get_L( Mesh& mesh);
    void get_L2( Mesh& mesh);
    bool particle_boundaries(double3& coord);
    void get_P();
    void get_Pr();
    void move_virt(Mesh& mesh,int timestep);
    void bound_resumption(const Particle& particle, const double3& r_new, const double3& p_new);

    void set_space_distribution();
    void set_pulse_distribution();
    void set_uniform_circle(int3 start, int3 end);
    void set_strict_uniform(int3 start, int3 end);
    void set_uniform(int3 start, int3 end);

    void add_uniform_cilinder(int numParts, double r0, double z0, double3 c);
protected:
    World &_world;
    double _mass;
    double _mpw; /*macroparticle weight*/
};



double PulseFromKev(double kev, double mass);

int get_num_of_type_particles(const std::vector<ParticlesArray> &Particles, const std::string& ParticlesType);
/// Ionization particles = electron(particles_e) +  ion (particles_i) 
void collision(const Mesh &mesh, const World& world ,std::vector<ParticlesArray> &Particles,int timestep);
//double getFieldEInLaser(double z, double r, const Region& domain, int timestep);
void reserve_Lmat(Mesh& mesh, ParticlesArray &sp);
#endif 
