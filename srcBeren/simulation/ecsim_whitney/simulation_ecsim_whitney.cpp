#include "simulation_ecsim_whitney.h"

#include <stdexcept>

#include "config/boundary_config.h"
#include "config/domain_config.h"
#include "particles/boris_pusher.h"
#include "runtime/profiler.h"
#include "whitney_shape.h"

using namespace ecsimWhitney;

void SimulationEcsimWhitney::init() {
    RECORD_TIMER;
    // Validate before allocating the mesh or creating output files. A
    // non-periodic axis is supported only when every face terminates the
    // orbit; reflecting walls need surface charge and are rejected below.
    Domain checked_domain = domain_from_config(system_config);
    BoundaryConditionHandler checked_bc = boundary_conditions_from_config(system_config, checked_domain);
    bool all_periodic = true;
    for (int a = 0; a < 3; ++a) {
        if (checked_domain.num_cells(a) < 4)
            throw std::invalid_argument("ecsim_whitney requires >=4 cells per axis");
        if (!checked_bc.is_periodic(a))
            all_periodic = false;
    }
    if (!all_periodic && !supports_wall_boundaries())
        throw std::invalid_argument("this Whitney predictor does not support non-periodic walls");
    if (all_periodic) {
        for (const auto& entry : system_config.at("Boundary_conditions"))
            if (entry.size() != 1 || !entry.contains("periodic"))
                throw std::invalid_argument("an all-periodic Whitney box requires periodic boundary entries");
    } else {
        for (const auto& entry : system_config.at("Boundary_conditions")) {
            const std::string type = entry.begin().key();
            if (type == "bphi" || type == "electron_reflection")
                throw std::invalid_argument("ecsim_whitney walls support absorbing open/ER0 faces only");
        }
        const auto& policy = checked_bc.path_policy();
        const bool cylinder = checked_domain.geom.use_cylinder;
        for (int a = 0; a < 3; ++a) {
            if (checked_bc.is_periodic(a))
                continue;
            // Cylinder geometry replaces the X/Y box faces by one absorbing
            // cylindrical boundary; the Z ends must still be absorbing.
            if (cylinder && a < 2) {
                if (!policy.clips(Face::CYLINDER))
                    throw std::invalid_argument("ecsim_whitney: cylinder side must have an absorbing boundary");
                continue;
            }
            const Face lower[3] = {Face::XMIN, Face::YMIN, Face::ZMIN};
            const Face upper[3] = {Face::XMAX, Face::YMAX, Face::ZMAX};
            if (!policy.clips(lower[a]) || !policy.clips(upper[a]))
                throw std::invalid_argument(
                    "ecsim_whitney: every face of a non-periodic axis must have an absorbing boundary");
        }
    }
    // CylinderDomain is supported: the cylindrical side is an absorbing wall.
    if (system_config.value("DampingType", "None") != "None")
        throw std::invalid_argument("ecsim_whitney requires no damping");
    const double dt = system_config.at("Dt").get<double>();
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("ecsim_whitney requires finite positive Dt");
    const std::string assembly = system_config.value("WhitneyAssembly", "reference");
    if (assembly != "reference" && assembly != "cell")
        throw std::invalid_argument("WhitneyAssembly must be reference or cell");
    parallel_point_ = assembly == "cell";
    const std::string kernel = system_config.value("WhitneyKernel", "linear");
    if (kernel == "linear")
        kernel_ = KernelOrder::Linear;
    else if (kernel == "quadratic")
        kernel_ = KernelOrder::Quadratic;
    else
        throw std::invalid_argument("WhitneyKernel must be linear or quadratic");
    SimulationEcsim::init();
    wall_boundaries_ = false;
    for (int a = 0; a < 3; ++a) {
        domain.periodic[a] = bc_handler.is_periodic(a);
        wall_boundaries_ = wall_boundaries_ || !domain.is_periodic(a);
    }
    sync_periodic(fieldE, domain);
    sync_periodic(fieldEn, domain);
    sync_periodic(fieldB, domain);
    sync_periodic(fieldBn, domain);
    sync_periodic(fieldBInit, domain);
    sync_periodic(fieldE_external, domain);
}

void SimulationEcsimWhitney::drift_particles_and_assemble_response() {
    RECORD_TIMER;
    fieldBFull = fieldB + fieldBInit;
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        // x^{n-1/2} -> x^{n+1/2}, with the ordinary ECSIM drift v^n.
        sp.move(dt_);
        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);
    }
    assemble_point_response();
}

void SimulationEcsimWhitney::assemble_point_response() {
    RECORD_TIMER;
    std::vector<Triplet> entries;
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        if (sp.is_neutral())
            continue;
        const double qm = sp.charge() / sp.mass(), qw = sp.charge() * sp.mpw();
        // Every point in a bucket uses the same 12 edge DOFs. Cells assemble
        // independently; thread-local currents avoid shared scatter writes.
        auto assemble_cell = [&](int cell, Field3d& current_grid, std::vector<Triplet>& matrix) {
            if (kernel_ == KernelOrder::Quadratic) {
                // The S2 stencil spans neighbouring cells, but within one cell
                // only eight index sets exist: the base on each axis is i-1 for
                // f<1/2 and i otherwise. Accumulate one dense block per base
                // instead of emitting 54x54 triplets per particle.
                static constexpr int slots = max_electric_weights * max_electric_weights;
                static thread_local std::vector<double> local;
                static thread_local std::array<ElectricStencil, 8> group_stencil;
                static thread_local std::array<bool, 8> used{};
                if (local.size() < 8 * slots)
                    local.resize(8 * slots);
                std::fill(local.begin(), local.begin() + 8 * slots, 0.0);
                used.fill(false);
                for (const auto& p : sp(cell)) {
                    const auto u = domain.to_cell_coordinates(p.coord);
                    int key = 0;
                    for (int a = 0; a < 3; ++a)
                        if (u[a] - std::floor(u[a]) < 0.5)
                            key |= 1 << a;
                    const auto stencil = electric_stencil(p.coord, domain, kernel_);
                    if (!used[key]) {
                        group_stencil[key] = stencil;
                        used[key] = true;
                    }
                    const auto b = (0.5 * dt_ * qm) * gather_magnetic(fieldBFull, p.coord, domain, kernel_);
                    const auto current = qw * magnetic_response(p.velocity, b);
                    for (int i = 0; i < stencil.size; ++i)
                        current_grid[stencil.weight[i].index] +=
                            stencil.weight[i].weight * current[stencil.weight[i].component];
                    Vector3R columns[3];
                    for (int c = 0; c < 3; ++c) {
                        Vector3R unit(0, 0, 0);
                        unit[c] = 1;
                        columns[c] = (0.25 * dt_ * dt_ * qw * qm) * magnetic_response(unit, b);
                    }
                    double* block = local.data() + key * slots;
                    for (int i = 0; i < stencil.size; ++i) {
                        const auto& row = stencil.weight[i];
                        for (int j = 0; j < stencil.size; ++j) {
                            const auto& col = stencil.weight[j];
                            block[i * max_electric_weights + j] +=
                                row.weight * columns[col.component][row.component] * col.weight;
                        }
                    }
                }
                for (int key = 0; key < 8; ++key) {
                    if (!used[key])
                        continue;
                    const auto& stencil = group_stencil[key];
                    const double* block = local.data() + key * slots;
                    for (int i = 0; i < stencil.size; ++i)
                        for (int j = 0; j < stencil.size; ++j) {
                            const double value = block[i * max_electric_weights + j];
                            if (value != 0)
                                matrix.emplace_back(stencil.weight[i].index, stencil.weight[j].index, value);
                        }
                }
                return;
            }
            double local[12][12] = {};
            std::array<EdgeWeight, 12> edges{};
            bool occupied = false;
            for (const auto& p : sp(cell)) {
                edges = electric_weights(p.coord, domain);
                occupied = true;
                const auto b = (0.5 * dt_ * qm) * gather_magnetic(fieldBFull, p.coord, domain, kernel_);
                const auto current = qw * magnetic_response(p.velocity, b);
                for (const auto& row : edges) current_grid[row.index] += row.weight * current[row.component];
                // L = dt^2 q^2 w/(4m) W^T alpha W, so
                // J_bar = Jp + (2/dt) L E_bar. Units match the ECSIM field solve.
                for (int c = 0; c < 3; ++c) {
                    Vector3R unit(0, 0, 0);
                    unit[c] = 1;
                    const auto column = (0.25 * dt_ * dt_ * qw * qm) * magnetic_response(unit, b);
                    for (int i = 0; i < 12; ++i)
                        for (int j = 4 * c; j < 4 * c + 4; ++j)
                            local[i][j] += edges[i].weight * column[edges[i].component] * edges[j].weight;
                }
            }
            if (occupied)
                for (int i = 0; i < 12; ++i)
                    for (int j = 0; j < 12; ++j)
                        if (local[i][j] != 0)
                            matrix.emplace_back(edges[i].index, edges[j].index, local[i][j]);
        };
        if (!parallel_point_) {
            for (int cell = 0; cell < sp.size(); ++cell) assemble_cell(cell, fieldJp, entries);
        } else {
            const int threads = omp_get_max_threads();
            std::vector<std::vector<Triplet>> matrices(threads);
            std::vector<Field3d> currents;
            for (int t = 0; t < threads; ++t) currents.push_back(Field3d::Zero(domain.size(), 3));
#pragma omp parallel num_threads(threads)
            {
                const int tid = omp_get_thread_num();
#pragma omp for schedule(static)
                for (int cell = 0; cell < sp.size(); ++cell) assemble_cell(cell, currents[tid], matrices[tid]);
            }
            size_t count = entries.size();
            for (const auto& matrix : matrices) count += matrix.size();
            entries.reserve(count);
            for (int t = 0; t < threads; ++t) {
                entries.insert(entries.end(), matrices[t].begin(), matrices[t].end());
                fieldJp += currents[t];
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
        // Mixed periodic/open boxes need both the alias copies and the
        // eliminated wall rows/columns; the handler applies them together.
        field_system.L = physical;
        bc_handler.apply_to_operator(field_system.L, domain);
    } else {
        // Maxwell stores periodic aliases as extra rows. Copy each physical
        // response row to its aliases, keeping columns on unique physical DOFs.
        // Applying the deposit-folding boundary operator here would double count.
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
        field_system.L.setFromTriplets(entries.begin(), entries.end());
    }
}

void SimulationEcsimWhitney::advance_particle_velocities() {
    RECORD_TIMER;
    const Field3d electric = fieldEp + fieldE_external;
    fieldBFull = fieldB + fieldBInit;
    for (auto& [name, ptr] : species) {
        auto& sp = *ptr;
        sp.currentOnGrid.setZero();
        if (sp.is_neutral())
            continue;
        const double qm = sp.charge() / sp.mass(), qw = sp.charge() * sp.mpw();
        auto push_cell = [&](int cell, Field3d& current_grid) {
            for (auto& p : sp(cell)) {
                const auto stencil = electric_stencil(p.coord, domain, kernel_);
                Vector3R gathered_electric(0, 0, 0);
                for (int i = 0; i < stencil.size; ++i) {
                    const auto& w = stencil.weight[i];
                    gathered_electric[w.component] += w.weight * electric[w.index];
                }
                const auto before = p.velocity;
                borisPusher::update_vEB(p, qm, gathered_electric, gather_magnetic(fieldBFull, p.coord, domain, kernel_),
                                        dt_);
                const auto current = (0.5 * qw) * (before + p.velocity);
                for (int i = 0; i < stencil.size; ++i) {
                    const auto& w = stencil.weight[i];
                    current_grid[w.index] += w.weight * current[w.component];
                }
            }
        };
        if (!parallel_point_) {
            for (int cell = 0; cell < sp.size(); ++cell) push_cell(cell, sp.currentOnGrid);
        } else {
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
    }
}
