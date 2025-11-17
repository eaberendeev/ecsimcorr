#ifndef PARTICLES_ARRAY_H_
#define PARTICLES_ARRAY_H_
#include <assert.h>

#include <functional>
#include <memory>

#include "Mesh.h"
#include "Particle.h"
#include "World.h"
#include "containers.h"
#include "random_generator.h"

typedef Eigen::Triplet<double> Trip;

class EnergySpectrum {
   public:
    EnergySpectrum(const std::vector<int>& spec, double minE,
                   double maxE)
        : minEnergy(minE), maxEnergy(maxE) {
        spectrum.resize(spec.size());
        for (size_t i = 0; i < spectrum.size(); ++i) {
            spectrum[i] = spec[i];
        }
    };
    EnergySpectrum(){};
    double minEnergy;
    double maxEnergy;
    std::vector<int> spectrum;
};

struct InterpolationWeights {
    alignas(32) int indices[2];
    alignas(32) double weights[2];
};

static inline InterpolationWeights compute_weights(double coord, double shift) {
    InterpolationWeights result;
    const double shifted_coord = coord - shift;
    const int index = static_cast<int>(shifted_coord);
    const double delta = shifted_coord - index;

    result.indices[0] = index;
    result.indices[1] = index + 1;
    result.weights[0] = 1.0 - delta;
    result.weights[1] = delta;

    return result;
}

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
    int3 cell;
    //#pragma omp declare simd linear(i : 1), notinbranch
    constexpr double& operator()(int i, int comp) {
        return shape[i + 2 * comp];
    }
};
double3 interpolateE_Chen(const Field3d& fieldE, const Node& node, ShapeK& sh,
                          ShapeK& sh_n);
double3 interpolateE(const Field3d& fieldE, const Node& node, ShapeK& no,
                     ShapeK& sh);
double3 interpolateB(const Field3d& fieldB, const Node& node, ShapeK& no,
                     ShapeK& sh);
double3 interpolateE_Chen(const Field3d& fieldE, ShapeK& sh, ShapeK& sh_n);
double3 interpolateB(const Field3d& fieldB, ShapeK& shape, ShapeK& shape05);
class ParticlesArray{

public:
 using ShapeFunction = double (*)(const double&);

 // We store p[articles by cells. In each cell we store vector of particles
 Array3D<std::vector<Particle>> particlesData;
 // how many particles in each cell. Need to be updated after move particles
 Array3D<int> countInCell;
 // struct ShapeK;
 void fill_shape(const Node& node, ShapeK& shape, bool shift) const;

 bool boundary_correction(double3& coord);
 bool boundary_correction(double3& coord, const int dim);
 Field3d densityOnGrid;
 Field3d currentOnGrid;

 std::vector<std::string> distSpace;
 std::vector<std::string> distPulse;
 std::string distType;

#ifdef SET_PARTICLE_IDS
    size_t ids;
#endif
    ParametersMap sortParameters;

    const double charge;
    const double density;
    double phasePXmin, phasePXmax;
    double kineticEnergy;
    double injectionEnergy;
    double lostEnergyZ;
    double lostEnergyXY;
    double lostParticlesZ;
    double lostParticlesXY;
    const std::string _name;
    const double3 temperature;
    const int NumPartPerCell;
    void delete_bounds();
    void add_particle(Particle &particle);
    void add_particles(std::vector<Particle>& particles);
    void save_init_coord();
    void save_init_velocity();
    void save_init_coord_and_velocity();
    void update_cells(const Domain& domain);
    void set_particles();

    void move_x(const double dt, Field3d& fieldJ);
    void move_y(const double dt, Field3d& fieldJ);

    void distribute_particles(const ParametersMap& parameters,
                              const Domain& domain, double timestep);
    std::vector<Particle> distribute_particles_in_space(
        const ParametersMap& parameters, ThreadRandomGenerator& randGenSpace);
    double distribute_particles_pulse(
        std::vector<Particle>& particles, const ParametersMap& parameters,
        ThreadRandomGenerator& randGenPulse);
    const std::string& name() const noexcept { return _name; }
    bool is_neutral() const noexcept { return charge == 0 || name() == "Neutrals"; }
    // delete particle if it lost the it cell

    void delete_particle_runtime(int ix, int iy, int iz, int ip) {
        countInCell(ix, iy, iz)--;
        int old_count = countInCell(ix, iy, iz);
        particlesData(ix, iy, iz)[ip] = particlesData(ix, iy, iz)[old_count];
        int lastParticle = particlesData(ix, iy, iz).size() - 1;
        if (old_count == lastParticle) {
            particlesData(ix, iy, iz).pop_back();
        } else {
            particlesData(ix, iy, iz)[old_count] =
                particlesData(ix, iy, iz)[lastParticle];
            particlesData(ix, iy, iz).pop_back();
        }
    }

    void delete_particle_runtime(int id_cell, int ip) {
        countInCell(id_cell)--;
        int old_count = countInCell(id_cell);
        particlesData(id_cell)[ip] = particlesData(id_cell)[old_count];
        int lastParticle = particlesData(id_cell).size() - 1;
        if (old_count == lastParticle) {
            particlesData(id_cell).pop_back();
        } else {
            particlesData(id_cell)[old_count] =
                particlesData(id_cell)[lastParticle];
            particlesData(id_cell).pop_back();
        }
    }

    void update_count_in_cell() {
        for (auto i = 0; i < countInCell.size().x(); i++) {
            for (auto j = 0; j < countInCell.size().y(); j++) {
                for (auto k = 0; k < countInCell.size().z(); k++) {
                    countInCell(i, j, k) = particlesData(i, j, k).size();
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
    virtual ~ParticlesArray() = default;

    void set_params_from_string(const std::string& line);
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
    double get_init_kinetic_energy() const;
    double3 get_kinetic_energy_component() const;
    double get_kinetic_energy(int dim) const;
    double get_kinetic_energy(int dim1, int dim2) const;
    double get_inject_energy() const { return injectionEnergy; }
    void glue_density_bound();
    void move(double dt);
    void prepare();

    void correctv_component(const Field3d& fieldE, const Field3d& fieldEp,
                            const Field3d& fieldEn, const Domain& domain,
                            const double dt);

    bool particle_boundaries(Particle& particle, const Domain& domain);
    void get_P();

    template <typename VelocityCalculator1, typename VelocityCalculator2>
    void calculate_pressure_component(Field3d& P,
                                      VelocityCalculator1 velocityCalc1, 
                                      VelocityCalculator2 velocityCalc2);

    EnergySpectrum calculate_energy_spectrum() const;
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
    // implicit methods
    void push_Chen(const Field3d& fieldE, const Field3d& fieldB, double dt);
    void calc_current_Chen(const double3& coord_start, const double3& coord_end,
                                           Field3d& fieldJ, const double dt);
    void updateJ_Chen(const double3 value, Field3d& fieldJ,
                                      const Node& node, ShapeK& sh,
                                      ShapeK& sh_n);
    double3 get_relative_coord(const double3& r) const {
        return double3(r.x() / xCellSize, r.y() / yCellSize, r.z() / zCellSize);
    }

    void fill_shape_from_coord(const double3& __r, ShapeK& shape, bool isShift) const;
    void fill_shape_from_voxel_and_coord(const int3& voxel, const double3& __r, ShapeK& shape,
                    bool isShift) const;
    void fill_shape(const int3& voxel, const double3& __r, ShapeK& shape) const;

    void updateJ_Chen(const double3 value, Field3d& fieldJ, ShapeK& sh, ShapeK& sh_n);



    void density_on_grid_update(ShapeType type = SHAPE);
    void predict_velocity(const Field3d& fieldE, const Field3d& fieldEp,
                          const Field3d& fieldB, const Domain& domain,
                          const double dt, ShapeType type = SHAPE);

    void move_and_calc_current(const double dt, Field3d& fieldJ,
                               ShapeType type = SHAPE);
    void predict_current(const Field3d& fieldB, Field3d& fieldJ,
                         const Domain& domain, const double dt,
                         ShapeType type = SHAPE);

    void fill_matrixL(Mesh& mesh, const Field3d& fieldB, const Domain& domain,
                      const double dt, ShapeType type = SHAPE);
    void fill_matrixL2(Mesh& mesh, const Field3d& fieldB, const Domain& domain,
                      const double dt, ShapeType type = SHAPE);

    double3 normalize_coord(const double3& coord) const {
        return double3(coord.x() / xCellSize, coord.y() / yCellSize,
                       coord.z() / zCellSize);
    }

   protected:

    template <ShapeFunction ShapeFn, int ShapeSize>
    void density_on_grid_update_impl();

    template <ShapeFunction ShapeFn, int ShapeSize>
    void move_and_calc_current_impl(const double dt, Field3d& fieldJ);

    void density_on_grid_update_impl_ngp();

    void fill_matrixL_impl_ngp(Mesh& mesh, const Field3d& fieldB,
                          const Domain& domain, const double dt);
    void fill_matrixL_impl_ngp2(Mesh& mesh, const Field3d& fieldB,
                                   const Domain& domain, const double dt);
    void fill_matrixL_impl_linear(Mesh& mesh, const Field3d& fieldB,
                          const Domain& domain, const double dt);
    void fill_matrixL_impl_linear2(Mesh& mesh, const Field3d& fieldB,
                                  const Domain& domain, const double dt);
    void predict_velocity_impl_ngp(const Field3d& fieldE, const Field3d& fieldEp,
                             const Field3d& fieldB, const Domain& domain,
                             const double dt);
    void predict_velocity_impl_linear(const Field3d& fieldE,
                                   const Field3d& fieldEp,
                                   const Field3d& fieldB, const Domain& domain,
                                   const double dt);
    void predict_current_impl_ngp(const Field3d& fieldB, Field3d& fieldJ,
                                  const Domain& domain, const double dt);
    void predict_current_impl_linear(const Field3d& fieldB, Field3d& fieldJ,
                                  const Domain& domain, const double dt);

    auto get_cell_index(const double3& coord) const {
        return std::array<int, 3>{int(coord.x() / xCellSize + GHOST_CELLS),
                                  int(coord.y() / yCellSize + GHOST_CELLS),
                                  int(coord.z() / zCellSize + GHOST_CELLS)};
    }
    double _mass;
    double _mpw; /*macroparticle weight*/
    double xCellSize;
    double yCellSize;
    double zCellSize;
    int xCellCount;
    int yCellCount;
    int zCellCount;
    Bounds bounds;
};

struct RadialVelocity {
    double operator()(const double3& coord, const double3& velocity,
                      double x0) const {
        double R = sqrt((coord.x() - x0) * (coord.x() - x0) +
                        (coord.y() - x0) * (coord.y() - x0));
        return ((coord.x() - x0) / R) * velocity.x() +
               ((coord.y() - x0) / R) * velocity.y();
    }
};

struct PhiVelocity {
    double operator()(const double3& coord, const double3& velocity,
                      double x0) const {
        double R = sqrt((coord.x() - x0) * (coord.x() - x0) +
                        (coord.y() - x0) * (coord.y() - x0));
        return -((coord.y() - x0) / R) * velocity.x() +
               ((coord.x() - x0) / R) * velocity.y();
    }
};

struct ZVelocity {
    double operator()([[maybe_unused]] const double3& coord,
                      const double3& velocity,
                      [[maybe_unused]] double x0) const {
        return velocity.z();
    }
};

typedef std::vector<std::unique_ptr<ParticlesArray>> Species;

double PulseFromKev(double kev, double mass);

int get_num_of_type_particles(const Species& Particles,
                              const std::string& ParticlesType);
/// Ionization particles = electron(particles_e) +  ion (particles_i)
void collision(const Mesh& mesh, Species& Particles, int timestep);
void reserve_Lmat(Mesh& mesh, ParticlesArray &sp);

#endif 
