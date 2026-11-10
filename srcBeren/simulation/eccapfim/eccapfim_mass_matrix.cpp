// Frozen-response mass matrix (dt/2) dJ/dE of the particles and its solve.
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "algorithms_eccapfim.h"
#include "runtime/profiler.h"
#include "simulation_eccapfim.h"

void SimulationEccapfim::add_particle_mass(const Vector3R& position, const Eigen::Matrix3d& response) {
    mass_points_[omp_get_thread_num()].push_back({position, response});
}

// Adds W T W^T of one particle. Rows lie on nodes cell-1 .. cell+2 of the
// particle cell, so particles of cells four apart never share a row.
void SimulationEccapfim::deposit_particle_mass(const Vector3R& position, const Eigen::Matrix3d& response) {
    // E_c at (i + 1/2, j, k) for c = x, and cyclic: linear weights along c on
    // half-integer nodes and quadratic weights across c on integer nodes.
    const Vector3R r = domain.to_cell_coordinates(position);
    const auto spline1 = [](double x) { return std::max(0.0, 1.0 - std::abs(x)); };
    const auto spline2 = [](double x) {
        x = std::abs(x);
        return x <= 0.5 ? 0.75 - x * x : (x < 1.5 ? 0.5 * (1.5 - x) * (1.5 - x) : 0.0);
    };

    struct Weight {
        int node[3];
        int row;
        double value;
    };

    Weight weights[3][18];
    int count[3] = {0, 0, 0};
    for (int c = 0; c < 3; ++c) {
        int first[3], width[3];
        double w[3][3];
        for (int a = 0; a < 3; ++a) {
            if (a == c) {
                first[a] = static_cast<int>(std::floor(r[a] - 0.5));
                width[a] = 2;
                for (int n = 0; n < 2; ++n) w[a][n] = spline1(r[a] - (first[a] + n + 0.5));
            } else {
                first[a] = static_cast<int>(std::floor(r[a] + 0.5)) - 1;
                width[a] = 3;
                for (int n = 0; n < 3; ++n) w[a][n] = spline2(r[a] - (first[a] + n));
            }
        }
        for (int i = 0; i < width[0]; ++i)
            for (int j = 0; j < width[1]; ++j)
                for (int k = 0; k < width[2]; ++k) {
                    const double value = w[0][i] * w[1][j] * w[2][k];
                    if (value == 0.0)
                        continue;
                    const int node[3] = {first[0] + i, first[1] + j, first[2] + k};
                    const int row = reduced_dof(node[0] + GHOST_CELLS, node[1] + GHOST_CELLS, node[2] + GHOST_CELLS, c);
                    weights[c][count[c]++] = {{node[0], node[1], node[2]}, row, value};
                }
    }
    for (int c1 = 0; c1 < 3; ++c1)
        for (int a = 0; a < count[c1]; ++a) {
            const Weight& wa = weights[c1][a];
            // Constrained rows keep the identity of the Maxwell operator.
            if (wa.row < 0 || electric_free_[wa.row] == 0)
                continue;
            for (int c2 = 0; c2 < 3; ++c2) {
                const double t = response(c1, c2) * wa.value;
                if (t == 0.0)
                    continue;
                for (int b = 0; b < count[c2]; ++b) {
                    const Weight& wb = weights[c2][b];
                    stencil_operator_.add_unsynchronized(wa.row, c1, c2, wb.node[0] - wa.node[0],
                                                         wb.node[1] - wa.node[1], wb.node[2] - wa.node[2],
                                                         t * wb.value);
                }
            }
        }
}

void SimulationEccapfim::assemble_mass_matrix(const Field3d& magnetic) {
    // Crank-Nicolson response: v_h = R(beta) (v^n + alpha E), alpha = dt q / 2m,
    // beta = alpha B, so (dt/2) dJ/dE = (dt/2)^2 q^2 w / m R(beta).
    for (const auto& [name, state] : kinetic_states_) {
        const auto& sp = *state.species;
        if (sp.is_neutral())
            continue;
        const double q = sp.charge(), m = sp.mass(), mpw = sp.mpw();
#pragma omp parallel for schedule(dynamic, 16)
        for (int cell = 0; cell < state.trial.capacity(); ++cell)
            for (size_t p = 0; p < state.trial(cell).size(); ++p) {
                // The orbit kernels average E along x^n -> x^{n+1}; its midpoint is the
                // closest point kernel. Orbits split by wall events are rare and use it too.
                const Vector3R midpoint = mass_matrix_point(state.initial(cell)[p].coord, state.trial(cell)[p].coord);
                const Vector3R B = algorithmsECCAPFIM::gather_magnetic_point(magnetic, midpoint, domain, &bc_handler);
                const Vector3R beta = B * (0.5 * dt_ * q / m);
                const double denom = 1.0 + beta.squaredNorm();
                const double a = 0.25 * dt_ * dt_ * q * q * mpw / m / denom;
                Eigen::Matrix3d T;
                T << a * (1 + beta.x() * beta.x()), a * (beta.z() + beta.x() * beta.y()),
                    a * (-beta.y() + beta.x() * beta.z()), a * (-beta.z() + beta.y() * beta.x()),
                    a * (1 + beta.y() * beta.y()), a * (beta.x() + beta.y() * beta.z()),
                    a * (beta.y() + beta.z() * beta.x()), a * (-beta.x() + beta.z() * beta.y()),
                    a * (1 + beta.z() * beta.z());
                add_particle_mass(algorithmsECCAPFIM::wrap_position(midpoint, domain, &bc_handler), T);
            }
    }
}

Vector3R SimulationEccapfim::mass_matrix_point(const Vector3R& start, const Vector3R& end) const {
    if (mass_matrix_at_end_)
        return end;
    // The trial end point may be wrapped through a periodic face: undo that shift.
    Vector3R displacement = end - start;
    for (int d = 0; d < 3; ++d)
        if (bc_handler.is_periodic(d)) {
            const double length = domain.num_cells(d) * domain.cell_size(d);
            displacement[d] -= length * std::round(displacement[d] / length);
        }
    return start + 0.5 * displacement;
}

// Deposits the collected particles without atomic updates: cells are
// coloured so that cells of one colour are at least four apart on every axis
// (also across periodic faces); a colour is processed in parallel.
void SimulationEccapfim::deposit_mass_points() {
    const Vector3I n = domain.num_cells();
    int colours[3], colour_of_first_remainder[3];
    for (int a = 0; a < 3; ++a) {
        const int full = n[a] >= 4 ? 4 * (n[a] / 4) : 0;
        colour_of_first_remainder[a] = full > 0 ? 4 : 0;
        colours[a] = colour_of_first_remainder[a] + (n[a] - full);
    }
    const auto colour = [&](int c, int a) {
        const int full = n[a] >= 4 ? 4 * (n[a] / 4) : 0;
        return c < full ? c % 4 : colour_of_first_remainder[a] + (c - full);
    };
    const int cells = n.x() * n.y() * n.z();
    const int total_colours = colours[0] * colours[1] * colours[2];
    // Counting sort of the points by colour, then by cell.
    std::vector<int> count(total_colours + 1, 0);
    std::vector<std::pair<int, const MassPoint*>> keyed;
    for (const auto& points : mass_points_)
        for (const auto& point : points) {
            const Vector3R r = domain.to_cell_coordinates(point.position);
            int c[3];
            for (int a = 0; a < 3; ++a) c[a] = std::clamp(static_cast<int>(std::floor(r[a])), 0, n[a] - 1);
            const int key = ((colour(c[0], 0) * colours[1] + colour(c[1], 1)) * colours[2] + colour(c[2], 2)) * cells +
                            (c[0] * n.y() + c[1]) * n.z() + c[2];
            keyed.emplace_back(key, &point);
        }
    std::sort(keyed.begin(), keyed.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    size_t begin = 0;
    while (begin < keyed.size()) {
        const int group = keyed[begin].first / cells;
        size_t end = begin;
        while (end < keyed.size() && keyed[end].first / cells == group) ++end;
        // Cell boundaries inside the colour group; each cell is one task.
        std::vector<size_t> starts{begin};
        for (size_t i = begin + 1; i < end; ++i)
            if (keyed[i].first != keyed[i - 1].first)
                starts.push_back(i);
        starts.push_back(end);
        const long tasks = static_cast<long>(starts.size()) - 1;
#pragma omp parallel for schedule(dynamic, 1)
        for (long task = 0; task < tasks; ++task)
            for (size_t i = starts[task]; i < starts[task + 1]; ++i)
                deposit_particle_mass(keyed[i].second->position, keyed[i].second->response);
        begin = end;
    }
}

// A = maxwell_ + (dt/2) M_p in the fixed particle pattern; GMRES with node
// block-Jacobi by default (xpic: PETSc GMRES, rtol 1e-6), or sparse LU.
void SimulationEccapfim::rebuild_mass_matrix(const Field3d& magnetic) {
    RECORD_TIMER;
    {
        RECORD_TIMER_NAMED("mass matrix assembly");
        stencil_operator_.reset();
        mass_points_.resize(omp_get_max_threads());
        for (auto& points : mass_points_) points.clear();
        assemble_mass_matrix(magnetic);
        deposit_mass_points();
    }
    timer::timer factorization_timer("mass matrix factorization");
    if (mass_matrix_lu_) {
        mass_preconditioner_ = stencil_operator_.to_sparse();
        mass_solver_.analyzePattern(mass_preconditioner_);
        mass_solver_.factorize(mass_preconditioner_);
        if (mass_solver_.info() != Eigen::Success)
            throw std::runtime_error("ECCAPFIM: mass-matrix factorization failed: " + mass_solver_.lastErrorMessage());
    } else {
        stencil_operator_.factor(mass_preconditioner_kind_);
    }
    mass_matrix_ready_ = true;
}

// Built once per step at the first use; the solvers clear the flag at every
// outer iteration when `mass_matrix_every_iteration` is set.
void SimulationEccapfim::ensure_mass_matrix(const Field3d& magnetic) {
    if (!mass_matrix_ready_)
        rebuild_mass_matrix(magnetic);
}

SimulationEccapfim::Vector SimulationEccapfim::solve_mass_system(const Vector& rhs) {
    RECORD_TIMER;
    Vector solution;
    if (mass_matrix_lu_) {
        solution = mass_solver_.solve(rhs);
        if (mass_solver_.info() != Eigen::Success)
            throw std::runtime_error("ECCAPFIM: mass-matrix solve failed");
    } else {
        int iterations = 0;
        if (!stencil_operator_.solve_gmres(rhs, solution, mass_matrix_tolerance_, 50, 2000, iterations))
            throw std::runtime_error("ECCAPFIM: mass-matrix GMRES did not converge");
        mass_solve_iterations_ += iterations;
    }
    if (!solution.allFinite())
        throw std::runtime_error("ECCAPFIM: nonfinite mass-matrix solution");
    return solution;
}

Eigen::SparseMatrix<double> SimulationEccapfim::mass_only_matrix() const {
    return stencil_operator_.to_sparse() - maxwell_;
}
