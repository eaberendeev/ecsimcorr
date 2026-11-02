#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "core/vector3.h"
#include "grid/geometry.h"

// Fixed for the lifetime of a species. Charge/mass and density use the code
// normalization; macro_weight = density / particles_per_cell (no cell volume).
// Restart restores particles, while these parameters come from the input.
struct SpeciesParameters {
    std::string name;
    double charge;
    double mass;
    double density;
    int particles_per_cell;
    double macro_weight;
};

class EnergySpectrum {
   public:
    EnergySpectrum(const std::vector<int>& spec, double minE, double maxE) : minEnergy(minE), maxEnergy(maxE) {
        spectrum.resize(spec.size());
        for (size_t i = 0; i < spectrum.size(); ++i) {
            spectrum[i] = spec[i];
        }
    }
    EnergySpectrum() {
    }
    double minEnergy;
    double maxEnergy;
    std::vector<int> spectrum;
};

struct BoundaryStepEvents {
    double lost_energy = 0;
    int lost_count = 0;
    int reflected_count = 0;
    double emitted_energy = 0;
    int emitted_count = 0;
};

// Reset before injection; measurements only read these events. Thread-local
// boundary events merge after the parallel classification phase. Damping belongs
// to the field step, not to a species. Recomputed rather than stored in recovery.
struct SpeciesStepState {
    double kinetic_energy_before_step = 0;
    double injection_energy = 0;
    int injection_count = 0;
    std::unordered_map<Face, BoundaryStepEvents> boundary;

    double total_lost_energy() const {
        double sum = 0;
        for (const std::pair<const Face, BoundaryStepEvents>& kv : boundary) sum += kv.second.lost_energy;
        return sum;
    }
    double total_lost_energy_z() const {
        double sum = 0;
        std::unordered_map<Face, BoundaryStepEvents>::const_iterator it = boundary.find(Face::ZMIN);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        it = boundary.find(Face::ZMAX);
        if (it != boundary.end())
            sum += it->second.lost_energy;
        return sum;
    }
    double total_lost_energy_xy() const {
        double sum = 0;
        std::unordered_map<Face, BoundaryStepEvents>::const_iterator it = boundary.find(Face::XMIN);
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
        std::unordered_map<Face, BoundaryStepEvents>::const_iterator it = boundary.find(Face::ZMIN);
        if (it != boundary.end())
            sum += it->second.lost_count;
        it = boundary.find(Face::ZMAX);
        if (it != boundary.end())
            sum += it->second.lost_count;
        return sum;
    }
    int total_lost_count_xy() const {
        int sum = 0;
        std::unordered_map<Face, BoundaryStepEvents>::const_iterator it = boundary.find(Face::XMIN);
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
        BoundaryStepEvents& s = boundary[face];
        s.lost_energy += energy;
        s.lost_count += 1;
    }
    void add_reflected(Face face) {
        boundary[face].reflected_count += 1;
    }
    void add_emitted(Face face, double energy) {
        BoundaryStepEvents& s = boundary[face];
        s.emitted_energy += energy;
        s.emitted_count += 1;
    }
    // Сливает диагностику из локального (поточного) буфера в общую: суммирует
    // поля инжекции и попарно все поля BoundaryStepEvents каждой грани.
    void merge_from(const SpeciesStepState& other) {
        injection_energy += other.injection_energy;
        injection_count += other.injection_count;
        for (const std::pair<const Face, BoundaryStepEvents>& kv : other.boundary) {
            BoundaryStepEvents& s = boundary[kv.first];
            s.lost_energy += kv.second.lost_energy;
            s.lost_count += kv.second.lost_count;
            s.reflected_count += kv.second.reflected_count;
            s.emitted_energy += kv.second.emitted_energy;
            s.emitted_count += kv.second.emitted_count;
        }
    }
    void begin_step(double kinetic_energy) {
        kinetic_energy_before_step = kinetic_energy;
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
