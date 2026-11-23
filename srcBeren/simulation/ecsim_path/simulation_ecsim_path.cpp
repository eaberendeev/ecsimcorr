#include "simulation_ecsim_path.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "diagnostics/Diagnostic.h"
#include "path_gather.h"
#include "path_push.h"
#include "runtime/profiler.h"

void SimulationEcsimPath::init_fields() {
    if (system_config.contains("EnergyCorrection"))
        throw std::invalid_argument("ecsim_path has no energy correction");
    SimulationEcsimCorr::init_fields();
    corrector_passes_ = system_config.value("CorrectorPasses", 1);
    if (corrector_passes_ < 1)
        throw std::invalid_argument("ecsim_path needs CorrectorPasses >= 1");
    magnetization_moment_.resize(domain.size(), 3);
    magnetization_current_.resize(domain.size(), 3);
    magnetization_current_.setZero();
}

void SimulationEcsimPath::prepare_step(int timestep) {
    SimulationEcsimCorr::prepare_step(timestep);
    if (!path_checked_) {
        warn_about_long_paths();
        path_checked_ = true;
    }
}

// Thermal path per step of the fastest charged species against the smallest cell.
void SimulationEcsimPath::warn_about_long_paths() const {
    const Vector3R h = domain.cell_size();
    const double cell = std::min(h.x(), std::min(h.y(), h.z()));
    for (const std::pair<const std::string, std::reference_wrapper<ParticlesArray>>& kv : charged_species) {
        const ParticlesArray& sp = kv.second.get();
        const double count = sp.get_total_num_of_particles();
        if (!(count > 0))
            continue;
        const double v2 = 2 * sp.get_kinetic_energy() / (sp.mass() * sp.mpw() * count);
        const double path = std::sqrt(v2 / 3) * dt_ / cell;
        if (path > 0.5)
            std::cerr << "ecsim_path: thermal path of " << sp.name() << " per step is " << path
                      << " cells (> 0.5): the scheme may be unstable, reduce Dt\n";
    }
}

void SimulationEcsimPath::predict_species_current(ParticlesArray& sp) {
    pathGather::predict_current<ShapeKernel<SHAPE_CH>::weight>(sp, fieldBFull, bc_handler.path_policy(sp, domain),
                                                               dt_, fieldJp);
}

void SimulationEcsimPath::after_first_half_step() {
    RECORD_TIMER;
    magnetization_moment_.setZero();
    const Field3d electric = fieldE + fieldE_external;
    for (const std::pair<const std::string, std::reference_wrapper<ParticlesArray>>& kv : charged_species)
        pathMagnetization::deposit(kv.second.get(), fieldBFull, electric, dt_, moments_[kv.first],
                                   magnetization_moment_);
    // Periodic aliases and halo deposits are folded into their owners and
    // copied back, as for a current (the fold does not depend on the layout).
    bc_handler.apply_to_fields(magnetization_moment_, FieldType::CURRENT, domain);
    spmv(maxwell.curlB, magnetization_moment_, magnetization_current_);
    // curl_B of the folded dM is complete on the owner edges: mask it like E
    // (a CURRENT fold would add the periodic aliases twice).
    bc_handler.apply_to_fields(magnetization_current_, FieldType::ELECTRIC, domain);
}

void SimulationEcsimPath::add_extra_current(Field3d& current) const {
    current += magnetization_current_;
}

void SimulationEcsimPath::push_species_and_deposit(ParticlesArray& sp, const Field3d& electric,
                                                   const Field3d& magnetic, SpeciesSecondHalfStep& result) {
    push_on_path_and_deposit<ShapeKernel<SHAPE_CH>::weight, ShapeKernel<SHAPE_CH>::deposit_half_width>(
        sp, electric, magnetic, dt_, sp.currentOnGrid, result.predictor_power,
        bc_handler.path_policy(sp, domain, &result.boundary_events), gather_passes_);
}

// The part a particle cannot give (v_perp too small) goes to a species factor.
void SimulationEcsimPath::after_magnetic_field_advance() {
    RECORD_TIMER;
    magnetization_work_ = 0;
    magnetization_remainder_ = 0;
    const Field3d change = fieldBn - fieldB;
    for (const std::pair<const std::string, std::reference_wrapper<ParticlesArray>>& kv : charged_species) {
        ParticlesArray& sp = kv.second.get();
        const pathMagnetization::Exchange result =
            pathMagnetization::exchange(sp, fieldBFull, change, dt_, moments_.at(kv.first));
        const double remainder = result.target - result.given;
        magnetization_work_ += result.target;
        magnetization_remainder_ += remainder;
        const double kinetic = sp.get_kinetic_energy();
        if (kinetic > 0 && remainder != 0) {
            const double lambda = std::sqrt(1 + remainder / kinetic);
            if (!std::isfinite(lambda) || lambda <= 0)
                throw std::runtime_error("invalid magnetization factor for " + sp.name());
#pragma omp parallel for schedule(dynamic, sp.omp_granularity())
            for (int cell = 0; cell < sp.size(); ++cell)
                for (Particle& p : sp(cell)) p.velocity = lambda * p.velocity;
        }
    }
}

// Each particle already gains the work of the field on its own path.
void SimulationEcsimPath::correct_species_energy(ParticlesArray&, double) {
}

void SimulationEcsimPath::diagnostic_values(Diagnostics& diagnostic) const {
    SimulationEcsimCorr::diagnostic_values(diagnostic);
    diagnostic.addValue("magnetizationWork", magnetization_work_);
    diagnostic.addValue("magnetizationRemainder", magnetization_remainder_);
}
