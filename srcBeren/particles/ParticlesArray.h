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

/**
 * @brief Storage for particle's coordinate - `r` (global, in PetscReal
 * units of dx, dy, dz), and cell - `g`, shifted cell g_s
 * (rounded, shifted by `shape_radius`).
 */
struct Node {
    double3 r;
    int3 g;
    int3 g05;

    Node(const double3& __r, const double3& cellSize) {
        set(__r, cellSize);
    }
    void set(const double3& __r, const double3& cellSize) {
        r = double3(__r.x() / cellSize.x(), __r.y() / cellSize.y(),
                    __r.z() / cellSize.z());

        g = int3(int(r.x() + 1) - 1, int(r.y() + 1) - 1, int(r.z() + 1) - 1);
        g05 = int3(int(r.x() + 0.5) - 1, int(r.y() + 0.5) - 1,
                   int(r.z() + 0.5) - 1);
    }
};
struct ShapeK {
    // const int dim = 3;
    alignas(64) double shape[SHAPE_SIZE * 3];
    //#pragma omp declare simd linear(i : 1), notinbranch
    constexpr double& operator()(int i, int comp) {
        return shape[i + 2 * comp];
    }
};

class ParticlesArray{

public:
// We store p[articles by cells. In each cell we store vector of particles
    Array3D< std::vector<Particle> > particlesData;
    //struct ShapeK;
    void fill_shape(const Node& node, ShapeK& shape,
                    bool shift) const;
    double3 interpolateE_Chen(const Field3d& fieldE, const Node& node,
                              ShapeK& sh, ShapeK& sh_n);
    double3 interpolateE(const Field3d& fieldE, const Node& node, ShapeK& no,
                         ShapeK& sh);
    double3 interpolateB(const Field3d& fieldB, const Node& node, ShapeK& no,
                        ShapeK& sh);
    void push_Chen(const Field3d& fieldE, const Field3d& fieldB,
                             double dt);
    bool boundary_correction(double3& coord);
    bool boundary_correction(double3& coord, const int dim);
    Array3D<double> densityOnGrid;
    Array2D<double> phaseOnGrid;
    Field3d currentOnGrid;
    Array3D<double> Pxx;
    Array3D<double> Pyy;
    Array3D<double> Pzz;
    std::vector<std::string> distSpace;
    std::vector<std::string> distPulse;
    std::string distType;

    const double charge;
    const double density;
    double phasePXmin, phasePXmax;
    double kineticEnergy;
    double injectionEnergy;
    double lostEnergy;
    const std::string _name;
    const double3 temperature;
    const int NumPartPerCell;
    void delete_bounds();
    void add_particle(const Particle &particle);
    void add_particles(const std::vector<Particle>& particles);
    void save_init_coord();
    void save_init_velocity();
    void save_init_coord_and_velocity();
    void calc_Esirkepov_current(const double dt, Field3d& fieldJ) const;
    void update_cells(const Domain& domain);
    void set_particles();
    void distribute_particles(const ParametersMap& parameters, double timestep);
    std::vector<Particle> distribute_particles_in_space(
        const ParametersMap& parameters, ThreadRandomGenerator& randGenSpace);
    double distribute_particles_pulse(
        std::vector<Particle>& particles, const ParametersMap& parameters,
        ThreadRandomGenerator& randGenPulse);
    const std::string& name() const noexcept { return _name; }

    // delete particle if it lost the it cell
    void delete_particle_runtime(int ix, int iy, int iz, int ip) {
        std::swap(particlesData(ix, iy, iz)[ip],
                  particlesData(ix, iy, iz).back());
        particlesData(ix, iy, iz).pop_back();
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

    double track_particle(double3& coord, double3& velocity,
                          const Field3d& fieldE, const Field3d& fieldB,
                          double dt, bool& intersect_bound);
    //bool in_exended_domain(const double3& coord);

    bool is_voxel_in_area(const int3& voxel);
    bool make_periodic_bound_force(double3& point);

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
