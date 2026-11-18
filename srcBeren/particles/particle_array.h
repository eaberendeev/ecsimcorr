#ifndef PARTICLES_ARRAY_H_
#define PARTICLES_ARRAY_H_
#include <assert.h>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

#include "Particle.h"
#include "grid/domain.h"
#include "include/config.h"
#include "nlohmann/json.hpp"
#include "species_state.h"
#include "utils/containers.h"
#include "utils/memory.h"

struct ParticlesStorage {
    ParticlesStorage(const Vector3I sizes)
        : sizes(sizes), storages((sizes.elements_product() + chunkSize - 1) / chunkSize) {
    }

    VectorView<Particle> operator()(int64_t ix) {
        Storage& storage = storages[ix / chunkSize];
        const int64_t ixInStorage = ix % chunkSize;
        return storage.getVectorView(ixInStorage);
    }

    VectorView<const Particle> operator()(int64_t ix) const {
        const Storage& storage = storages[ix / chunkSize];
        const int64_t ixInStorage = ix % chunkSize;
        return storage.getVectorView(ixInStorage);
    }

    VectorView<Particle> operator()(int64_t i, int64_t j, int64_t k) {
        assert(i >= 0 && i < sizes[0]);
        assert(j >= 0 && j < sizes[1]);
        assert(k >= 0 && k < sizes[2]);
        const int64_t ix = getCellIx(i, j, k);
        return this->operator()(ix);
    }

    VectorView<const Particle> operator()(int64_t i, int64_t j, int64_t k) const {
        assert(i >= 0 && i < sizes[0]);
        assert(j >= 0 && j < sizes[1]);
        assert(k >= 0 && k < sizes[2]);
        const int64_t ix = getCellIx(i, j, k);
        return this->operator()(ix);
    }

    void addParticle(int64_t i, int64_t j, int64_t k, const Particle& particle) {
        const int ix = getCellIx(i, j, k);
        Storage& storage = storages[ix / chunkSize];
        storage.addParticle(ix % chunkSize, particle);
    }

    void removeParticle(int64_t i, int64_t j, int64_t k, int64_t particleLocalIx) {
        const int ix = getCellIx(i, j, k);
        Storage& storage = storages[ix / chunkSize];
        storage.removeParticle(ix % chunkSize, particleLocalIx);
    }

    int64_t capacity() const {
        return sizes.elements_product();
    }

    Vector3I size() const {
        return sizes;
    }

   private:
    int64_t getCellIx(int64_t i, int64_t j, int64_t k) const {
        return i * sizes[1] * sizes[2] + j * sizes[2] + k;
    }

    struct Storage {
        Storage() : particles(0), outerIndexes(chunkSize + 1), allocatedSize(0), usedSize(0) {
            std::fill_n(outerIndexes.get(), chunkSize + 1, 0);
        }

        std::vector<Particle> particles;
        SmartPtr<int> outerIndexes;
        int64_t allocatedSize;
        int64_t usedSize;

        void addParticle(int64_t ix, const Particle& particle) {
            particles.insert(particles.begin() + outerIndexes[ix], particle);
            for (int i = ix + 1; i < chunkSize + 1; ++i) {
                outerIndexes[i] += 1;
            }
        }

        void removeParticle(int64_t ix, int64_t particleLocalIx) {
            particles.erase(particles.begin() + outerIndexes[ix] + particleLocalIx);
            for (int i = ix + 1; i < chunkSize + 1; ++i) {
                outerIndexes[i] -= 1;
            }
        }

        VectorView<Particle> getVectorView(int64_t ix) {
            const int64_t start = outerIndexes[ix];
            const int64_t end = outerIndexes[ix + 1];
            if (end == start) {
                return VectorView<Particle>(nullptr, 0);
            }
            VectorView<Particle> res(&particles[start], end - start);
            return res;
        }

        VectorView<const Particle> getVectorView(int64_t ix) const {
            const int64_t start = outerIndexes[ix];
            const int64_t end = outerIndexes[ix + 1];
            if (end == start) {
                return VectorView<const Particle>(nullptr, 0);
            }
            VectorView<const Particle> res(&particles[start], end - start);
            return res;
        }
    };

    static constexpr int64_t chunkSize = 32;

    Vector3I sizes;
    std::vector<Storage> storages;
};

class ParticlesArray {
   private:
    const SpeciesParameters parameters_;
    // One shared simulation domain; its address and grid layout outlive this species.
    const Domain& domain_;
    ParticlesStorage cells_;

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

    ParticlesStorage& cells() {
        return cells_;
    }

    const ParticlesStorage& cells() const {
        return cells_;
    }

    // void swap_and_pop_particle(int ix, int iy, int iz, int index) {
    //     swap_and_pop_particle(domain_.sind(ix, iy, iz), index);
    // }

    // void swap_and_pop_particle(int cell, int index) {
    //     std::vector<Particle>& particles = cells_(cell);
    //     std::swap(particles[index], particles.back());
    //     particles.pop_back();
    // }

    virtual ~ParticlesArray() = default;

    double mass() const {
        return parameters_.mass;
    }

    double mpw() const {
        return parameters_.macro_weight;
    }

    VectorView<Particle> operator()(int i) {
        return cells_(i);
    }

    const VectorView<const Particle> operator()(int i) const {
        return cells_(i);
    }

    VectorView<Particle> operator()(int i, int j, int k) {
        return cells_(i, j, k);
    }

    const VectorView<const Particle> operator()(int i, int j, int k) const {
        return cells_(i, j, k);
    }

    Vector3I sizes() const {
        return cells_.size();
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
