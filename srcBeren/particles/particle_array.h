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
#include "species_state.h"

class ParticlesArray {
   private:
    const SpeciesParameters parameters_;
    // One shared simulation domain; its address and grid layout outlive this species.
    const Domain& domain_;
    Array3D<std::vector<Particle>> cells_;

   public:
    ParticlesArray(const nlohmann::json& config, const Domain& domain);
    ParticlesArray(const nlohmann::json& config, const Domain&& domain) = delete;

    ParticlesArray(ParticlesArray&& other) = default;
    ParticlesArray(const ParticlesArray& other) = delete;

    // Accumulated current and energy/events belong to the numerical step.
    // begin_step() snapshots kinetic energy and resets current/events before injection;
    // reading diagnostics does not consume them. These are not restart payloads.
    Field3d currentOnGrid;

    SpeciesStepState step_state;
    void add_particle(const Particle& particle);
    void add_particles(const std::vector<Particle>& particles);
    // Drift and boundary handling must already have placed particles in their
    // current or neighboring cell. Insertion/removal invalidate boundary snapshots.
    void update_cells(const Domain& domain);

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
        return cells_;
    }
    const Array3D<std::vector<Particle>>& cells() const {
        return cells_;
    }

    void swap_and_pop_particle(int ix, int iy, int iz, int index) {
        swap_and_pop_particle(domain_.sind(ix, iy, iz), index);
    }
    void swap_and_pop_particle(int cell, int index) {
        std::vector<Particle>& particles = cells_(cell);
        std::swap(particles[index], particles.back());
        particles.pop_back();
    }

    virtual ~ParticlesArray() = default;

    double mass() const {
        return parameters_.mass;
    }

    double mpw() const {
        return parameters_.macro_weight;
    }
    std::vector<Particle>& operator()(int i) {
        return cells_(i);
    }

    const std::vector<Particle>& operator()(int i) const {
        return cells_(i);
    }
    std::vector<Particle>& operator()(int i, int j, int k) {
        return cells_(i, j, k);
    }
    const std::vector<Particle>& operator()(int i, int j, int k) const {
        return cells_(i, j, k);
    }
    int size() const {
        return cells_.capacity();
    }
    ParticleMoments diagnostic_moments() const;
    double get_kinetic_energy() const;
    Vector3R get_kinetic_energy_component() const;
    double get_kinetic_energy(int dim) const;
    double get_kinetic_energy(int dim1, int dim2) const;
    double get_inject_energy() const {
        return step_state.injection_energy;
    }
    void move(double dt);
    void begin_step();

    EnergySpectrum calculate_energy_spectrum() const;
    int get_total_num_of_particles() const {
        int count = 0;
        for (int cell = 0; cell < size(); ++cell) count += cells_(cell).size();
        return count;
    }

    const Domain& get_domain() const {
        return domain_;
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
