#ifndef PARTICLES_ARRAY_H_
#define PARTICLES_ARRAY_H_
#include <assert.h>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

#include "Particle.h"
#include "config.h"
#include "containers.h"
#include "grid/domain.h"
#include "nlohmann/json.hpp"
#include "particle_storage.h"
#include "particles_distribution_collection.h"
#include "random.h"
#include "species_parameters.h"

class EnergySpectrum {
   public:
    EnergySpectrum(const std::vector<int>& spec, double minE, double maxE) : minEnergy(minE), maxEnergy(maxE) {
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
   private:
    const SpeciesParameters parameters_;
    ParticleStorage storage_;

   public:
    ParticlesArray(const nlohmann::json& config, const Domain& domain);
    ParticlesArray(const nlohmann::json& config, const Domain&& domain) = delete;

    ParticlesArray(ParticlesArray&& other) = default;
    ParticlesArray(const ParticlesArray& other) = delete;

    // Accumulated current and energy/events belong to the numerical step.
    // prepare() snapshots kinetic energy and resets current/events before injection;
    // reading diagnostics does not consume them. These are not restart payloads.
    Field3d currentOnGrid;

    double kineticEnergy = 0;
    SpeciesDiagStats diag;
    void add_particle(const Particle& particle) {
        storage_.add_particle(particle);
    }
    void add_particles(const std::vector<Particle>& particles) {
        storage_.add_particles(particles);
    }
    void update_cells(const Domain& domain) {
        storage_.update_cells(domain);
    }

    const std::string& name() const noexcept {
        return parameters_.name;
    }
    bool is_neutral() const noexcept {
        return charge() == 0 || name() == "Neutrals";
    }

    const SpeciesParameters& parameters() const noexcept {
        return parameters_;
    }
    double charge() const noexcept {
        return parameters_.charge;
    }
    double density() const noexcept {
        return parameters_.density;
    }
    int particles_per_cell() const noexcept {
        return parameters_.particles_per_cell;
    }
    Array3D<std::vector<Particle>>& cells() {
        return storage_.cells();
    }
    const Array3D<std::vector<Particle>>& cells() const {
        return storage_.cells();
    }

    void swap_and_pop_particle(int ix, int iy, int iz, int index) {
        storage_.swap_and_pop_particle(get_domain().sind(ix, iy, iz), index);
    }
    void swap_and_pop_particle(int id_cell, int index) {
        storage_.swap_and_pop_particle(id_cell, index);
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

    double mass() const {
        return parameters_.mass;
    }

    double mpw() const {
        return parameters_.macro_weight;
    }
    std::vector<Particle>& operator()(int i) {
        return cells()(i);
    }

    const std::vector<Particle>& operator()(int i) const {
        return cells()(i);
    }
    int size() const {
        return storage_.size();
    }
    ParticleMoments diagnostic_moments() const;
    double get_kinetic_energy() const;
    Vector3R get_kinetic_energy_component() const;
    double get_kinetic_energy(int dim) const;
    double get_kinetic_energy(int dim1, int dim2) const;
    double get_inject_energy() const {
        return diag.injection_energy;
    }
    void move(double dt);
    void prepare();

    EnergySpectrum calculate_energy_spectrum() const;
    int get_total_num_of_particles() const {
        return storage_.particle_count();
    }

    const Domain& get_domain() const {
        return storage_.domain();
    }

    int omp_granularity(const int nthr = omp_get_max_threads()) const {
        return std::min(512, std::max(size() / nthr, 1));
    }
};

// Keyed by ParticlesArray::name() (species name from config)
using Species = std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>;

// Returns nullptr if not found
ParticlesArray* find_species(Species& species, const std::string& name);
const ParticlesArray* find_species(const Species& species, const std::string& name);

#endif
