#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "drift_kinetic_push.h"
#include "grid/domain.h"
#include "utils/containers.h"

namespace driftKinetic {

// Cell-bucketed guiding-center storage. The cell key is the same storage bin
// used by the kinetic ParticlesArray; periodic wrapping during a push is the
// caller's responsibility, while this container performs a full rebin when
// positions change by more than one cell.
class DriftKineticParticles {
   public:
    DriftKineticParticles(const Domain& domain, std::string name, double charge, double mass, double weight = 1.0)
        : domain_(domain), name_(std::move(name)), charge_(charge), mass_(mass), weight_(weight),
          cells_(domain.size()) {
        if (!(mass_ > 0.0) || !std::isfinite(mass_))
            throw std::invalid_argument("drift-kinetic species requires positive finite mass");
        if (!std::isfinite(charge_))
            throw std::invalid_argument("drift-kinetic species requires finite charge");
    }

    void add(const GuidingCenter& particle) {
        const Vector3I cell = domain_.storage_cell_index(particle.position);
        if (!domain_.grid.contains_cell_index(cell))
            throw std::out_of_range("drift-kinetic particle is outside the stored grid");
        cells_(cell.x(), cell.y(), cell.z()).push_back(particle);
    }

    void update_cells() {
        Array3D<std::vector<GuidingCenter>> cells(domain_.size());
        for (int cell = 0; cell < cells_.capacity(); ++cell)
            for (const auto& particle : cells_(cell)) {
                const Vector3I target = domain_.storage_cell_index(particle.position);
                if (!domain_.grid.contains_cell_index(target))
                    throw std::out_of_range("drift-kinetic particle is outside the stored grid");
                cells(target.x(), target.y(), target.z()).push_back(particle);
            }
        cells_ = std::move(cells);
    }

    const std::string& name() const {
        return name_;
    }
    double charge() const {
        return charge_;
    }
    double mass() const {
        return mass_;
    }
    double weight() const {
        return weight_;
    }
    bool is_neutral() const {
        return charge_ == 0.0;
    }
    const Domain& domain() const {
        return domain_;
    }

    int size() const {
        return cells_.capacity();
    }
    std::vector<GuidingCenter>& operator()(int cell) {
        return cells_(cell);
    }
    const std::vector<GuidingCenter>& operator()(int cell) const {
        return cells_(cell);
    }
    std::vector<GuidingCenter>& operator()(int i, int j, int k) {
        return cells_(i, j, k);
    }
    const std::vector<GuidingCenter>& operator()(int i, int j, int k) const {
        return cells_(i, j, k);
    }

    int get_total_num_of_particles() const {
        int total = 0;
        for (int cell = 0; cell < cells_.capacity(); ++cell) total += static_cast<int>(cells_(cell).size());
        return total;
    }

    std::unique_ptr<DriftKineticParticles> clone() const {
        auto copy = std::make_unique<DriftKineticParticles>(domain_, name_, charge_, mass_, weight_);
        copy->cells_ = cells_;
        return copy;
    }

    double parallel_kinetic_energy() const {
        double energy = 0.0;
        for (int cell = 0; cell < cells_.capacity(); ++cell)
            for (const auto& particle : cells_(cell))
                energy += 0.5 * mass_ * particle.v_parallel * particle.v_parallel;
        return energy;
    }

   private:
    const Domain& domain_;
    std::string name_;
    double charge_ = 0.0;
    double mass_ = 1.0;
    double weight_ = 1.0;
    Array3D<std::vector<GuidingCenter>> cells_;
};

// Binary recovery file of one species: a header (format tag, version,
// marker count) followed by the GuidingCenter records. Reading appends to
// `particles` and fails on a missing, foreign or truncated file.
void write_guiding_centers(const DriftKineticParticles& particles, const std::string& filename);
void read_guiding_centers(DriftKineticParticles& particles, const std::string& filename);

}   // namespace driftKinetic
