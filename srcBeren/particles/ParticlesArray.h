#ifndef PARTICLES_ARRAY_H_
#define PARTICLES_ARRAY_H_
#include <assert.h>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

#include "Mesh.h"
#include "Particle.h"
#include "boris_pusher.h"
#include "config.h"
#include "containers.h"
#include "decompose_esirkepov_current.h"
#include "grid/domain.h"
#include "nlohmann/json.hpp"
#include "particles_distribution_collection.h"
#include "random.h"
#include "row_block.h"

typedef Eigen::Triplet<double> Trip;

class EnergySpectrum {
   public:
    EnergySpectrum(const std::vector<int>& spec, double minE, double maxE) : minEnergy(minE), maxEnergy(maxE) {
        spectrum.resize(spec.size());
        for (size_t i = 0; i < spectrum.size(); ++i) {
            spectrum[i] = spec[i];
        }
    };
    EnergySpectrum() {};
    double minEnergy;
    double maxEnergy;
    std::vector<int> spectrum;
};

class ParticlesArray;

struct PerFaceStats {
    double lost_energy = 0;
    int lost_count = 0;
    int reflected_count = 0;
    double emitted_energy = 0;
    int emitted_count = 0;
};

struct SpeciesDiagStats {
    double injection_energy = 0;
    int injection_count = 0;
    std::unordered_map<Face, PerFaceStats> boundary;
    bool is_cylinder_geometry = false;

    double total_lost_energy() const {
        double sum = 0;
        for (const auto& kv : boundary) sum += kv.second.lost_energy;
        return sum;
    }
    double total_lost_energy_z() const {
        double sum = 0;
        auto it = boundary.find(Face::ZMIN);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        it = boundary.find(Face::ZMAX);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        return sum;
    }
    double total_lost_energy_xy() const {
        double sum = 0;
        auto it = boundary.find(Face::XMIN);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        it = boundary.find(Face::XMAX);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        it = boundary.find(Face::YMIN);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        it = boundary.find(Face::YMAX);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        it = boundary.find(Face::CYLINDER);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        return sum;
    }
    int total_lost_count_z() const {
        int sum = 0;
        auto it = boundary.find(Face::ZMIN);
        if (it != boundary.end())
            sum += it->second.lost_count;
        it = boundary.find(Face::ZMAX);
        if (it != boundary.end())
            sum += it->second.lost_count;
        return sum;
    }
    int total_lost_count_xy() const {
        int sum = 0;
        auto it = boundary.find(Face::XMIN);
        if (it != boundary.end())
            sum += it->second.lost_count;
        it = boundary.find(Face::XMAX);
        if (it != boundary.end())
            sum += it->second.lost_count;
        it = boundary.find(Face::YMIN);
        if (it != boundary.end())
            sum += it->second.lost_count;
        it = boundary.find(Face::YMAX);
        if (it != boundary.end())
            sum += it->second.lost_count;
        it = boundary.find(Face::CYLINDER);
        if (it != boundary.end())
            sum += it->second.lost_count;
        return sum;
    }
    void add_loss(Face face, double energy) {
        auto& s = boundary[face];
        s.lost_energy += energy;
        s.lost_count += 1;
    }
    void add_reflected(Face face) {
        boundary[face].reflected_count += 1;
    }
    void add_emitted(Face face, double energy) {
        auto& s = boundary[face];
        s.emitted_energy += energy;
        s.emitted_count += 1;
    }
    // Сливает диагностику из локального (поточного) буфера в общую: суммирует
    // поля инжекции и попарно все поля PerFaceStats каждой грани.
    void merge_from(const SpeciesDiagStats& other) {
        injection_energy += other.injection_energy;
        injection_count += other.injection_count;
        for (const auto& kv : other.boundary) {
            auto& s = boundary[kv.first];
            s.lost_energy += kv.second.lost_energy;
            s.lost_count += kv.second.lost_count;
            s.reflected_count += kv.second.reflected_count;
            s.emitted_energy += kv.second.emitted_energy;
            s.emitted_count += kv.second.emitted_count;
        }
    }
    void clear() {
        injection_energy = 0;
        injection_count = 0;
        boundary.clear();
    }
};

struct ParticleMoments {
    Vector3R kinetic{0, 0, 0};
    int64_t count = 0;
    double energy() const {
        return kinetic.x() + kinetic.y() + kinetic.z();
    }
};

class ParticlesArray {
   public:
    using ShapeFunction = double (*)(const double&);

    ParticlesArray(const nlohmann::json& config, const Domain& domain);
    ParticlesArray(const nlohmann::json& config, const Domain&& domain) = delete;

    ParticlesArray(ParticlesArray&& other) = default;
    ParticlesArray(const ParticlesArray& other) = delete;

    // We store p[articles by cells. In each cell we store vector of particles
    Array3D<std::vector<Particle>> particlesData;

    Field3d currentOnGrid;

    const double charge;
    const double density;
    double phasePXmin, phasePXmax;
    double kineticEnergy = 0;
    SpeciesDiagStats diag;
    const std::string name_;
    int NumPartPerCell;
    void add_particle(const Particle& particle);
    void add_particles(const std::vector<Particle>& particles);
    void update_cells(const Domain& domain);
    void set_particles();

    const std::string& name() const noexcept {
        return name_;
    }
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

    double distribute_initial_particles(const std::vector<std::unique_ptr<IDistribution>>& distributions,
                                        const Domain& domain, bool co_locate_species = true);
    double inject_particles_step(std::vector<std::unique_ptr<IDistribution>>& distributions, int timestep,
                                 const Domain& domain, double dt, bool co_locate_species = true);
    double add_particles_from_distribution(IDistribution& dist, LehmerEngine& rng_space, LehmerEngine& rng_momentum,
                                           const Domain& domain, double dt, bool check_boundaries);
    void add_distribution(const nlohmann::json& config, const std::string& type);
    const std::vector<std::unique_ptr<IDistribution>>& get_initial_distributions() const {
        return initialDistributions_;
    }
    const std::vector<std::unique_ptr<IDistribution>>& get_injection_distributions() const {
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
    double mass() const {
        return mass_;
    }

    double mpw() const {
        return mpw_;
    }
    std::vector<Particle>& operator()(int i) {
        return particlesData(i);
    }

    const std::vector<Particle>& operator()(int i) const {
        return particlesData(i);
    }
    int size() const {
        return particlesData.size().x() * particlesData.size().y() * particlesData.size().z();
    }
    ParticleMoments diagnostic_moments() const;
    double get_kinetic_energy() const;
    Vector3R get_kinetic_energy_component() const;
    double get_kinetic_energy(int dim) const;
    double get_kinetic_energy(int dim1, int dim2) const;
    double get_inject_energy() const {
        return diag.injection_energy;
    }
    void glue_density_bound();
    void move(double dt);
    void prepare();

    template <typename VelocityCalculator1, typename VelocityCalculator2>
    void calculate_pressure_component(Field3d& P, VelocityCalculator1 velocityCalc1,
                                      VelocityCalculator2 velocityCalc2) const;

    EnergySpectrum calculate_energy_spectrum() const;
    int get_total_num_of_particles() const {
        int count = 0;
        for (auto k = 0; k < size(); ++k) {
            count += particlesData(k).size();
        }

        return count;
    }
    template <ShapeType ST = SHAPE>
    void density_on_grid_update(Field3d& density) const;
    template <ShapeType ST = SHAPE>
    void density_on_grid_update_reference(Field3d& density) const;

    template <ShapeType ST = SHAPE>
    void fill_matrixL2_Reference(Mesh& mesh, const Field3d& fieldB, const Domain& domain, const double dt) const;
    template <ShapeType ST>
    void fill_matrixL2(Mesh& mesh, const Field3d& fieldB, const Domain& domain, const double dt,
                       SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>& rowBlocksGlobal) const;
    const Domain& get_domain() const {
        return domain_;
    }

    int omp_granularity(const int nthr = omp_get_max_threads()) const {
        return std::min(512, std::max(size() / nthr, 1));
    }

   protected:
    template <ShapeFunction ShapeFn, int ShapeSize>
    void density_on_grid_update_impl(Field3d& density) const;
    template <ShapeFunction ShapeFn, int ShapeSize>
    void density_on_grid_update_impl_reference(Field3d& density) const;

    void fill_matrixL_impl_linear2_Reference(Mesh& mesh, const Field3d& fieldB, const Domain& domain,
                                             const double dt) const;
    void fill_matrixL_impl_linear2(const Mesh& mesh, const Field3d& fieldB, const Domain& domain, const double dt,
                                   SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>& rowBlocksThrLocal) const;

    double mass_;
    double mpw_; /*macroparticle weight*/
    // All species use the same simulation geometry and mesh, without independent copies.
    // Borrowed from Simulation/test fixture; the domain outlives this species
    // and its grid layout stays fixed while particle cells exist.
    const Domain& domain_;
    nlohmann::json config_;
};

// The cylindrical basis is undefined exactly on the axis; choose the x/y
// basis there so moments remain finite and the pressure trace is preserved.
inline Vector3R cylindrical_velocity(const Vector3R& coord, const Vector3R& velocity, const Vector3R& center) {
    const double x = coord.x() - center.x(), y = coord.y() - center.y();
    const double radius2 = x * x + y * y;
    if (radius2 == 0)
        return velocity;
    const double inverse_radius = 1.0 / std::sqrt(radius2);
    return {(x * velocity.x() + y * velocity.y()) * inverse_radius,
            (-y * velocity.x() + x * velocity.y()) * inverse_radius, velocity.z()};
}

inline double radial_velocity(const Vector3R& coord, const Vector3R& velocity, const Vector3R& center) {
    return cylindrical_velocity(coord, velocity, center).x();
}
inline double phi_velocity(const Vector3R& coord, const Vector3R& velocity, const Vector3R& center) {
    return cylindrical_velocity(coord, velocity, center).y();
}
inline double z_velocity(const Vector3R&, const Vector3R& velocity, const Vector3R&) {
    return velocity.z();
}

// Keyed by ParticlesArray::name() (species name from config)
using Species = std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>;

double PulseFromKev(double kev, double mass);

// Returns nullptr if not found
ParticlesArray* find_species(Species& species, const std::string& name);
const ParticlesArray* find_species(const Species& species, const std::string& name);

/// Ionization particles = electron(particles_e) +  ion (particles_i)
void collision(const Mesh& mesh, Species& Particles, int timestep);
void reserve_Lmat(Mesh& mesh, ParticlesArray& sp);

#endif
