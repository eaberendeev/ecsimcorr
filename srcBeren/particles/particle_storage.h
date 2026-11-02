#pragma once

#include "Particle.h"
#include "containers.h"
#include "grid/domain.h"

// Owns the existing AoS cell vectors; no species parameters, sources or grid
// moments. The borrowed domain must outlive this storage and keep its layout.
// Insertion/removal/redistribution may invalidate particle addresses. Boundary
// snapshots must be consumed before any such operation.
class ParticleStorage {
   public:
    explicit ParticleStorage(const Domain& domain) : cells_(domain.size()), domain_(domain) {
    }
    ParticleStorage(const Domain&& domain) = delete;
    ParticleStorage(ParticleStorage&&) = default;
    ParticleStorage(const ParticleStorage&) = delete;

    Array3D<std::vector<Particle>>& cells() {
        return cells_;
    }
    const Array3D<std::vector<Particle>>& cells() const {
        return cells_;
    }
    const Domain& domain() const {
        return domain_;
    }
    int size() const {
        return cells_.capacity();
    }

    void add_particle(const Particle& particle);
    void add_particles(const std::vector<Particle>& particles);
    void swap_and_pop_particle(int cell, int index) {
        auto& particles = cells_(cell);
        std::swap(particles[index], particles.back());
        particles.pop_back();
    }
    int particle_count() const {
        int count = 0;
        for (int cell = 0; cell < size(); ++cell) count += cells_(cell).size();
        return count;
    }
    // Coordinates already include drift and physical boundary handling.
    // Every particle must stay in its current cell or move to a neighbor.
    void update_cells(const Domain& domain);

   private:
    Array3D<std::vector<Particle>> cells_;
    const Domain& domain_;
};
