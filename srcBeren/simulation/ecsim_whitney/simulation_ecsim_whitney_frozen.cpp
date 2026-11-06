#include "simulation_ecsim_whitney_frozen.h"

#include <atomic>
#include <exception>

using namespace ecsimWhitney;

void SimulationEcsimWhitneyFrozen::init() {
    const auto config = system_config.value("WhitneyOrbit", nlohmann::json::object());
    if (!config.is_object())
        throw std::invalid_argument("WhitneyOrbit must be an object");
    for (auto it = config.begin(); it != config.end(); ++it)
        if (it.key() != "predictor" && it.key() != "assembly")
            throw std::invalid_argument("unknown WhitneyOrbit option: " + it.key());
    const std::string predictor = config.value("predictor", "ballistic");
    if (predictor == "ballistic")
        predictor_ = OrbitPredictor::Ballistic;
    else if (predictor == "old_fields")
        predictor_ = OrbitPredictor::OldFields;
    else if (predictor == "half_drift")
        predictor_ = OrbitPredictor::HalfDrift;
    else if (predictor == "centered_fields")
        predictor_ = OrbitPredictor::CenteredFields;
    else
        throw std::invalid_argument(
            "WhitneyOrbit.predictor must be ballistic, old_fields, half_drift or centered_fields");
    const std::string assembly = config.value("assembly", "reference");
    if (assembly != "reference" && assembly != "cell")
        throw std::invalid_argument("WhitneyOrbit.assembly must be reference or cell");
    cell_assembly_ = assembly == "cell";
    SimulationEcsimWhitney::init();
}

void SimulationEcsimWhitneyFrozen::append_response_matrix(std::vector<Triplet>& entries,
                                                          const std::vector<ParticleResponse>& responses, size_t begin,
                                                          size_t end, double qm, double qw,
                                                          std::vector<int>& local_index) const {
    RECORD_TIMER_NAMED("frozen cell response matrix");
    // The union of orbit edge DOFs need not fit the point-shape 12x12 block.
    // Fold duplicate contributions within a cell before Eigen's global sort;
    // at large particles/cell this avoids hundreds of identical triplets.
    std::vector<int> dofs;
    if (cell_assembly_) {
        for (size_t p = begin; p < end; ++p)
            for (const auto& w : responses[p].electric)
                if (local_index[w.index] < 0) {
                    local_index[w.index] = static_cast<int>(dofs.size());
                    dofs.push_back(w.index);
                }
    }
    // Long predicted orbits retain the sparse reference path instead of
    // allocating an unbounded dense cell matrix.
    const size_t n = dofs.size() <= 256 ? dofs.size() : 0;
    std::vector<double> local(n * n, 0);
    for (size_t p = begin; p < end; ++p) {
        const auto& response = responses[p];
        const auto& weights = response.electric;
        const auto b = (0.5 * dt_ * qm) * response.magnetic;
        // Product of averaged gathers: L ~ Wbar^T alpha Wbar.
        // Averaging W^T alpha W instead would not reproduce the Boris push.
        Vector3R columns[3];
        for (int c = 0; c < 3; ++c) {
            Vector3R unit(0, 0, 0);
            unit[c] = 1;
            columns[c] = (0.25 * dt_ * dt_ * qw * qm) * magnetic_response(unit, b);
        }
        if (n == 0) {
            // Simple per-particle reference, including all off-diagonal
            // magnetic-response components and periodic owner indices.
            for (const auto& row : weights)
                for (const auto& col : weights) {
                    const double value = row.weight * columns[col.component][row.component] * col.weight;
                    if (value != 0)
                        entries.emplace_back(row.index, col.index, value);
                }
        } else {
            for (const auto& row : weights) {
                auto* local_row = local.data() + local_index[row.index] * n;
                const double alpha_row[3] = {row.weight * columns[0][row.component],
                                             row.weight * columns[1][row.component],
                                             row.weight * columns[2][row.component]};
                for (const auto& col : weights)
                    local_row[local_index[col.index]] += alpha_row[col.component] * col.weight;
            }
        }
    }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (local[i * n + j] != 0)
                entries.emplace_back(dofs[i], dofs[j], local[i * n + j]);
    for (int dof : dofs) local_index[dof] = -1;
}

void SimulationEcsimWhitneyFrozen::drift_particles_and_assemble_response() {
    RECORD_TIMER;
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        // Keep the actual ECSIM drift unchanged for a controlled comparison.
        sp.move(dt_);
        // With walls, keep absorbed charged particles until after the velocity
        // push: their truncated orbit still contributes to this step's J/L, and
        // the response order must match the particle buckets. Neutral species
        // have no response, so they can be boundary-processed immediately.
        if (!wall_boundaries_ || sp.is_neutral())
            bc_handler.apply_to_particles(sp, species, domain);
        sp.update_cells(domain);
    }
    assemble_frozen_response();
}

void SimulationEcsimWhitneyFrozen::assemble_frozen_response() {
    response_ready_ = false;
    std::erase_if(response_, [&](const auto& entry) {
        const auto it = species.find(entry.first);
        return it == species.end() || it->second->is_neutral();
    });
    fieldBFull = fieldB + fieldBInit;
    const Field3d old_electric = fieldE + fieldE_external;
    std::vector<Triplet> entries;
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        if (sp.is_neutral())
            continue;
        auto& responses = response_[name];
        const double qm = sp.charge / sp.mass(), qw = sp.charge * sp.mpw();
        const PathPolicy path_policy = wall_boundaries_ ? bc_handler.path_policy(sp, domain) : PathPolicy::none();
        // Shared physical kernel for the serial reference and independent
        // cell workers. Each worker supplies its own current deposit buffer.
        auto particle_response = [&](const Particle& p, Field3d& current_grid, ParticleResponse& response,
                                     OrbitWeights& weights, std::vector<double>& cuts) {
            // Reconstruct a periodic lift of the pre-drift coordinate. Never
            // join wrapped endpoints by a chord across the whole box.
            const auto start = p.coord - dt_ * p.velocity;
            const auto free_path =
                predict_path(start, p.velocity, qm, dt_, old_electric, fieldBFull, domain, predictor_);
            ParticleFateResult fate;
            const FrozenPath path =
                wall_boundaries_ ? clip_ballistic_path(free_path, p, path_policy, domain, fate) : free_path;
            if (cell_assembly_)
                orbit_weights_cell(path, domain, weights, cuts);
            else
                weights = orbit_weights(path, domain);
            response.removed = fate.fate == ParticleFate::Removed;
            const auto magnetic = gather_weights(fieldBFull, weights.magnetic);
            const auto b = (0.5 * dt_ * qm) * magnetic;
            const auto current = qw * magnetic_response(p.velocity, b);
            for (const auto& row : weights.electric) current_grid[row.index] += row.weight * current[row.component];
            // Only the electric gather and Bbar are needed after assembly.
            response.electric.assign(weights.electric.begin(), weights.electric.end());
            response.magnetic = magnetic;
        };
        if (!cell_assembly_) {
            responses.clear();
            std::vector<int> unused_index;
            OrbitWeights weights;
            std::vector<double> cuts;
            for (int cell = 0; cell < sp.size(); ++cell) {
                const size_t begin = responses.size();
                for (const auto& p : sp.particlesData(cell)) {
                    responses.emplace_back();
                    particle_response(p, fieldJp, responses.back(), weights, cuts);
                }
                append_response_matrix(entries, responses, begin, responses.size(), qm, qw, unused_index);
            }
        } else {
            // Prefix ranges preserve particle traversal order without shared
            // vector push_back or pointers invalidated by reallocation.
            std::vector<size_t> offsets(sp.size() + 1, 0);
            for (int cell = 0; cell < sp.size(); ++cell)
                offsets[cell + 1] = offsets[cell] + sp.particlesData(cell).size();
            responses.resize(offsets.back());
            const int threads = omp_get_max_threads();
            std::vector<std::vector<Triplet>> thread_entries(threads);
            std::vector<Field3d> thread_current;
            for (int t = 0; t < threads; ++t) thread_current.push_back(Field3d::Zero(domain.size(), 3));
            std::exception_ptr failure;
            std::atomic<bool> failed{false};
#pragma omp parallel num_threads(threads)
            {
                const int tid = omp_get_thread_num();
                auto& current = thread_current[tid];
                auto& matrix = thread_entries[tid];
                OrbitWeights weights;
                std::vector<double> cuts;
                // Only touched DOFs are reset between cells; there is no
                // assumed coordinate stencil for the union of orbit shapes.
                std::vector<int> local_index(3 * domain.total_size(), -1);
#pragma omp for schedule(static)
                for (int cell = 0; cell < sp.size(); ++cell) {
                    if (failed.load(std::memory_order_relaxed))
                        continue;
                    try {
                        size_t index = offsets[cell];
                        {
                            RECORD_TIMER_NAMED("frozen cell orbit weights and current");
                            for (const auto& p : sp.particlesData(cell))
                                particle_response(p, current, responses[index++], weights, cuts);
                        }
                        append_response_matrix(matrix, responses, offsets[cell], offsets[cell + 1], qm, qw,
                                               local_index);
                    } catch (...) {
#pragma omp critical(whitney_response_failure)
                        {
                            if (!failure)
                                failure = std::current_exception();
                        }
                        failed.store(true, std::memory_order_relaxed);
                    }
                }
            }
            if (failure)
                std::rethrow_exception(failure);
            size_t count = entries.size();
            for (const auto& matrix : thread_entries) count += matrix.size();
            entries.reserve(count);
            for (int t = 0; t < threads; ++t) {
                entries.insert(entries.end(), thread_entries[t].begin(), thread_entries[t].end());
                fieldJp += thread_current[t];
            }
        }
    }
    if (wall_boundaries_)
        bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);
    else
        sync_periodic(fieldJp, domain);
    filtered_current_power_ = 0;
    const int dofs = 3 * domain.total_size();
    Operator physical(dofs, dofs);
    physical.setFromTriplets(entries.begin(), entries.end());
    entries.clear();
    if (wall_boundaries_) {
        // The boundary handler applies periodic aliases and open/absorbing
        // row/column constraints together, including mixed periodic/wall boxes.
        mesh.Lmat2 = physical;
        bc_handler.apply_to_operator(mesh.Lmat2, domain);
    } else {
        const auto size = domain.size();
        for (int i = 0; i < size[0]; ++i)
            for (int j = 0; j < size[1]; ++j)
                for (int k = 0; k < size[2]; ++k)
                    for (int c = 0; c < 3; ++c) {
                        const int row = domain.vind(i, j, k, c);
                        const int source = domain.vind(periodic_index(i, 0, domain), periodic_index(j, 1, domain),
                                                       periodic_index(k, 2, domain), c);
                        for (Operator::InnerIterator it(physical, source); it; ++it)
                            entries.emplace_back(row, it.col(), it.value());
                    }
        mesh.Lmat2.setFromTriplets(entries.begin(), entries.end());
    }
    response_ready_ = true;
}

void SimulationEcsimWhitneyFrozen::advance_particle_velocities() {
    RECORD_TIMER;
    if (!response_ready_)
        throw std::runtime_error("frozen Whitney velocity push requires a newly assembled response");
    const Field3d electric = fieldEp + fieldE_external;
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        sp.currentOnGrid.setZero();
        if (sp.is_neutral())
            continue;
        const auto& responses = response_.at(name);
        std::vector<size_t> offsets(sp.size() + 1, 0);
        for (int cell = 0; cell < sp.size(); ++cell) offsets[cell + 1] = offsets[cell] + sp.particlesData(cell).size();
        if (offsets.back() != responses.size())
            throw std::runtime_error("frozen Whitney particle ordering changed before velocity push");
        const double qm = sp.charge / sp.mass(), qw = sp.charge * sp.mpw();
        auto push_cell = [&](int cell, Field3d& current_grid) {
            size_t index = offsets[cell];
            for (auto& p : sp.particlesData(cell)) {
                const auto& response = responses[index++];
                if (response.removed) {
                    // The absorbed orbit already contributes Jp + (2/dt)L E to
                    // Maxwell. Update the velocity at the wall but do not deposit
                    // a second current: the loss energy must include this work.
                    borisPusher::update_vEB(p, qm, gather_weights(electric, response.electric), response.magnetic, dt_);
                    continue;
                }
                const auto before = p.velocity;
                borisPusher::update_vEB(p, qm, gather_weights(electric, response.electric), response.magnetic, dt_);
                const auto current = (0.5 * qw) * (before + p.velocity);
                for (const auto& w : response.electric) current_grid[w.index] += w.weight * current[w.component];
            }
        };
        if (!cell_assembly_) {
            for (int cell = 0; cell < sp.size(); ++cell) push_cell(cell, sp.currentOnGrid);
        } else {
            // The response is immutable after the solve. Each worker owns its
            // particle range and deposits the adjoint work current locally.
            const int threads = omp_get_max_threads();
            std::vector<Field3d> currents;
            for (int t = 0; t < threads; ++t) currents.push_back(Field3d::Zero(domain.size(), 3));
#pragma omp parallel num_threads(threads)
            {
                const int tid = omp_get_thread_num();
#pragma omp for schedule(static)
                for (int cell = 0; cell < sp.size(); ++cell) push_cell(cell, currents[tid]);
            }
            for (const auto& current : currents) sp.currentOnGrid += current;
        }
        if (wall_boundaries_)
            bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
        else
            sync_periodic(sp.currentOnGrid, domain);
        if (wall_boundaries_)
            bc_handler.apply_to_particles(sp, species, domain);
    }
    response_ready_ = false;
}
