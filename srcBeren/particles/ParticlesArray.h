#ifndef PARTICLES_ARRAY_H_
#define PARTICLES_ARRAY_H_
#include <assert.h>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

#include "Mesh.h"
#include "Particle.h"
#include "World.h"
#include "boris_pusher.h"
#include "containers.h"
#include "decompose_esirkepov_current.h"
#include "nlohmann/json.hpp"
#include "particles_distribution_collection.h"
#include "random_generator.h"
#include "config.h"

typedef Eigen::Triplet<double> Trip;

class EnergySpectrum {
   public:
    EnergySpectrum(const std::vector<int>& spec, double minE, double maxE)
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
    Vector3R r;
    Vector3I g;
    Vector3I g05;

    Node(const Vector3R& __r, const Vector3R& cellSize) { set(__r, cellSize); }
    void set(const Vector3R& __r, const Vector3R& cellSize) {
        r = Vector3R(__r.x() / cellSize.x(), __r.y() / cellSize.y(),
                     __r.z() / cellSize.z());

        g = Vector3I(int(r.x() + 1) - 1, int(r.y() + 1) - 1,
                     int(r.z() + 1) - 1);
        g05 = Vector3I(int(r.x() + 0.5) - 1, int(r.y() + 0.5) - 1,
                       int(r.z() + 0.5) - 1);
    }
};
struct ShapeK {
    // const int dim = 3;
    alignas(64) double shape[SHAPE_SIZE * 3];
    Vector3I cell;
    // #pragma omp declare simd linear(i : 1), notinbranch
    constexpr double& operator()(int i, int comp) {
        return shape[i + 2 * comp];
    }
};
Vector3R interpolateE_Chen(const Field3d& fieldE, const Node& node, ShapeK& sh,
                           ShapeK& sh_n);
Vector3R interpolateE(const Field3d& fieldE, const Node& node, ShapeK& no,
                      ShapeK& sh);
Vector3R interpolateB(const Field3d& fieldB, const Node& node, ShapeK& no,
                      ShapeK& sh);
Vector3R interpolateE_Chen(const Field3d& fieldE, ShapeK& sh, ShapeK& sh_n);
Vector3R interpolateB(const Field3d& fieldB, ShapeK& shape, ShapeK& shape05);

class ParticlesArray {
   public:
    using ShapeFunction = double (*)(const double&);

    ParticlesArray(const nlohmann::json& config, const Domain& domain);

    // We store p[articles by cells. In each cell we store vector of particles
    Array3D<std::vector<Particle>> particlesData;

    // struct ShapeK;
    void fill_shape(const Node& node, ShapeK& shape, bool shift) const;

    Field3d densityOnGrid;
    Field3d currentOnGrid;

    const double charge;
    const double density;
    double phasePXmin, phasePXmax;
    double kineticEnergy;
    double injectionEnergy;
    double lostEnergyZ;
    double lostEnergyXY;
    double lostParticlesZ;
    double lostParticlesXY;
    const std::string name_;
    int NumPartPerCell;
    void add_particle(const Particle& particle);
    void add_particles(const std::vector<Particle>& particles);
    void save_init_coord();
    void save_init_velocity();
    void save_init_coord_and_velocity();
    void update_cells(const Domain& domain);
    void set_particles();

    void move_x(const double dt, Field3d& fieldJ);
    void move_y(const double dt, Field3d& fieldJ);

    const std::string& name() const noexcept { return name_; }
    bool is_neutral() const noexcept {
        return charge == 0 || name() == "Neutrals";
    }

    void swap_and_pop_particle(int ix, int iy, int iz, int index) {
        auto& cell = particlesData(ix, iy, iz);
        std::swap(cell[index], cell.back());
        cell.pop_back();
    }

    void swap_and_pop_particle(int id_cell, int index) {
        auto& cell = particlesData(id_cell);
        std::swap(cell[index], cell.back());
        cell.pop_back();
    }

    std::vector<std::unique_ptr<IDistribution>> initialDistributions_;
    std::vector<std::unique_ptr<IDistribution>> injectionDistributions_;

    void initialize_distributions(const nlohmann::json& config);
    void make_second_emission(const Particle& particle);

    double distribute_initial_particles(
        const std::vector<std::unique_ptr<IDistribution>>& distributions,
        const Domain& domain);
    double inject_particles_step(
        std::vector<std::unique_ptr<IDistribution>>& distributions,
        int timestep, const Domain& domain, double dt);
    double add_particles_from_distribution(IDistribution& dist,
                                           ThreadRandomGenerator& rng_space,
                                           ThreadRandomGenerator& rng_momentum,
                                           const Domain& domain, double dt,
                                           bool check_boundaries);
    void add_distribution(const nlohmann::json& config,
                          const std::string& type);
    const std::vector<std::unique_ptr<IDistribution>>&
    get_initial_distributions() const {
        return initialDistributions_;
    }
    const std::vector<std::unique_ptr<IDistribution>>&
    get_injection_distributions() const {
        return injectionDistributions_;
    }
    std::vector<std::unique_ptr<IDistribution>>& get_initial_distributions() {
        return initialDistributions_;
    }
    std::vector<std::unique_ptr<IDistribution>>& get_injection_distributions() {
        return injectionDistributions_;
    }
    virtual ~ParticlesArray() = default;

    void phase_on_grid_update(const Domain& domain);
    void inject(int timestep);
    void update(Mesh& mesh, int timestep);
    double mass() const { return mass_; }

    double mpw() const { return mpw_; }
    std::vector<Particle>& operator()(int i) { return particlesData(i); }

    const std::vector<Particle>& operator()(int i) const {
        return particlesData(i);
    }
    int size() const {
        return particlesData.size().x() * particlesData.size().y() *
               particlesData.size().z();
    }
    double get_kinetic_energy() const;
    double get_init_kinetic_energy() const;
    Vector3R get_kinetic_energy_component() const;
    double get_kinetic_energy(int dim) const;
    double get_kinetic_energy(int dim1, int dim2) const;
    double get_inject_energy() const { return injectionEnergy; }
    void glue_density_bound();
    void move(double dt);
    void prepare();

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
    double track_particle(Vector3R& coord, Vector3R& velocity,
                          const Field3d& fieldE, const Field3d& fieldB,
                          double dt, bool& intersect_bound);
    // bool in_exended_domain(const Vector3R& coord);

    bool is_voxel_in_area(const Vector3I& voxel);
    bool make_periodic_bound_force(Vector3R& point);
    // implicit methods
    void push_Chen(const Field3d& fieldE, const Field3d& fieldB, double dt);
    void calc_current_Chen(const Vector3R& coord_start,
                           const Vector3R& coord_end, Field3d& fieldJ,
                           const double dt);
    void updateJ_Chen(const Vector3R value, Field3d& fieldJ, const Node& node,
                      ShapeK& sh, ShapeK& sh_n);

    void fill_shape_from_coord(const Vector3R& __r, ShapeK& shape,
                               bool isShift) const;
    void fill_shape_from_voxel_and_coord(const Vector3I& voxel,
                                         const Vector3R& __r, ShapeK& shape,
                                         bool isShift) const;
    void fill_shape(const Vector3I& voxel, const Vector3R& __r,
                    ShapeK& shape) const;

    void updateJ_Chen(const Vector3R value, Field3d& fieldJ, ShapeK& sh,
                      ShapeK& sh_n);

    void density_on_grid_update(ShapeType type = SHAPE);

    void move_and_calc_current(const double dt, Field3d& fieldJ,
                               ShapeType type = SHAPE);

    void fill_matrixL(Mesh& mesh, const Field3d& fieldB, const Domain& domain,
                      const double dt, ShapeType type = SHAPE);
    void fill_matrixL2(Mesh& mesh, const Field3d& fieldB, const Domain& domain,
                       const double dt, ShapeType type = SHAPE);
    const auto& get_domain() const { return domain_; }

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

    double mass_;
    double mpw_; /*macroparticle weight*/
    Domain domain_;
    nlohmann::json config_;
};

struct RadialVelocity {
    double operator()(const Vector3R& coord, const Vector3R& velocity,
                      double x0) const {
        double R = sqrt((coord.x() - x0) * (coord.x() - x0) +
                        (coord.y() - x0) * (coord.y() - x0));
        return ((coord.x() - x0) / R) * velocity.x() +
               ((coord.y() - x0) / R) * velocity.y();
    }
};

struct PhiVelocity {
    double operator()(const Vector3R& coord, const Vector3R& velocity,
                      double x0) const {
        double R = sqrt((coord.x() - x0) * (coord.x() - x0) +
                        (coord.y() - x0) * (coord.y() - x0));
        return -((coord.y() - x0) / R) * velocity.x() +
               ((coord.x() - x0) / R) * velocity.y();
    }
};

struct ZVelocity {
    double operator()([[maybe_unused]] const Vector3R& coord,
                      const Vector3R& velocity,
                      [[maybe_unused]] double x0) const {
        return velocity.z();
    }
};

// Keyed by ParticlesArray::name() (species name from config)
using Species =
    std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>;

double PulseFromKev(double kev, double mass);

// Returns nullptr if not found
ParticlesArray* find_species(Species& species, const std::string& name);
const ParticlesArray* find_species(const Species& species,
                                   const std::string& name);

/// Ionization particles = electron(particles_e) +  ion (particles_i)
void collision(const Mesh& mesh, Species& Particles, int timestep);
void reserve_Lmat(Mesh& mesh, ParticlesArray& sp);

#endif
