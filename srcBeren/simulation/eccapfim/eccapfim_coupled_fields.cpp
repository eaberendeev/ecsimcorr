// Coupled E+B nonlinear solve: Maxwell residual, mass-matrix preconditioner,
// Anderson, NGMRES and Newton-Krylov iterations.
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

#include "algorithms_eccapfim.h"
#include "runtime/profiler.h"
#include "simulation_eccapfim.h"

SimulationEccapfim::Vector SimulationEccapfim::pack6(const Field3d& electric, const Field3d& magnetic) const {
    const auto n = unknown_range_.end - unknown_range_.start;
    const auto first = unknown_range_.start;
    Vector result(6 * n.elements_product());
    int p = 0;
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k) {
                for (int c = 0; c < 3; ++c) result[p++] = electric(i + first.x(), j + first.y(), k + first.z(), c);
                for (int c = 0; c < 3; ++c) result[p++] = magnetic(i + first.x(), j + first.y(), k + first.z(), c);
            }
    return result;
}

void SimulationEccapfim::unpack6(const Vector& vector, Field3d& electric, Field3d& magnetic) const {
    const auto n = unknown_range_.end - unknown_range_.start;
    const auto first = unknown_range_.start;
    int p = 0;
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k) {
                for (int c = 0; c < 3; ++c) electric(i + first.x(), j + first.y(), k + first.z(), c) = vector[p++];
                for (int c = 0; c < 3; ++c) magnetic(i + first.x(), j + first.y(), k + first.z(), c) = vector[p++];
            }
    algorithmsECCAPFIM::copy_periodic_nodes(electric, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(magnetic, domain, &bc_handler);
    bc_handler.apply_to_fields(electric, FieldType::ELECTRIC, domain);
    bc_handler.apply_to_fields(magnetic, FieldType::MAGNETIC, domain);
}

SimulationEccapfim::Vector SimulationEccapfim::coupled_residual(const Vector& state) {
    if (!state.allFinite())
        throw std::runtime_error("ECCAPFIM: nonfinite coupled field iterate");
    unpack6(state, fieldEp, fieldBHalf);
    const Field3d full_electric = fieldEp + fieldE_external;
    const Field3d full_magnetic = fieldBHalf + fieldBInit;
    ++residual_evaluations_;
    evaluate_current(full_electric, full_magnetic);

    // F_E = E_h - E^n + dt/2 (J - curl_B B_h),
    // F_B = B_h - B^n + dt/2 curl_E E_h.
    // `fieldJ` already contains the kinetic current and curl(M) for the
    // drift-kinetic backend.
    Field3d residual_e = fieldEp - fieldE + (0.5 * dt_) * (fieldJ - (maxwell.curlB * fieldBHalf));
    Field3d residual_b = fieldBHalf - fieldB + (0.5 * dt_) * (maxwell.curlE * fieldEp);
    algorithmsECCAPFIM::copy_periodic_nodes(residual_e, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(residual_b, domain, &bc_handler);
    bc_handler.apply_to_fields(residual_e, FieldType::ELECTRIC, domain);
    bc_handler.apply_to_fields(residual_b, FieldType::MAGNETIC, domain);
    return pack6(residual_e, residual_b);
}

SimulationEccapfim::Vector SimulationEccapfim::jacobian_vector(const Vector& state, const Vector& direction,
                                                               const Vector& base_residual) {
    const double direction_norm = direction.norm();
    if (direction_norm == 0.0)
        return Vector::Zero(base_residual.size());
    const double epsilon = 1e-7 * (1.0 + state.norm()) / direction_norm;
    const Vector perturbed = state + epsilon * direction;
    return (coupled_residual(perturbed) - base_residual) / epsilon;
}

SimulationEccapfim::Vector SimulationEccapfim::apply_right_preconditioner(const Vector& vector) {
    if (!use_mass_matrix_)
        return vector;
    // fieldBHalf holds the last evaluated iterate: B^n at the first call.
    ensure_mass_matrix(fieldBHalf + fieldBInit);
    Field3d field_e(domain.size(), 3), field_b(domain.size(), 3);
    unpack6(vector, field_e, field_b);
    const Field3d rhs_e = field_e + (0.5 * dt_) * (maxwell.curlB * field_b);
    const Vector rhs_vector = pack(rhs_e);
    const Vector correction_e = solve_mass_system(rhs_vector);
    Field3d correction_e_field(domain.size(), 3);
    unpack(correction_e, correction_e_field);
    Field3d correction_b = field_b - (0.5 * dt_) * (maxwell.curlE * correction_e_field);
    algorithmsECCAPFIM::copy_periodic_nodes(correction_e_field, domain, &bc_handler);
    algorithmsECCAPFIM::copy_periodic_nodes(correction_b, domain, &bc_handler);
    bc_handler.apply_to_fields(correction_e_field, FieldType::ELECTRIC, domain);
    bc_handler.apply_to_fields(correction_b, FieldType::MAGNETIC, domain);
    return pack6(correction_e_field, correction_b);
}

SimulationEccapfim::Vector SimulationEccapfim::gmres(const Vector& state, const Vector& rhs, int& iterations,
                                                     bool& ok) {
    const int n = rhs.size();
    iterations = 0;
    ok = true;
    if (n == 0)
        return Vector::Zero(0);
    const double beta = rhs.norm();
    if (beta == 0.0)
        return Vector::Zero(n);

    const int max_iterations = std::min(coupled_gmres_max_, n);
    std::vector<Vector> basis(max_iterations + 1);
    Eigen::MatrixXd hessenberg = Eigen::MatrixXd::Zero(max_iterations + 1, max_iterations);
    basis[0] = rhs / beta;
    const Vector base_residual = -rhs;   // rhs = -F(state).
    Vector solution = Vector::Zero(n);
    double linear_residual = beta;

    for (int j = 0; j < max_iterations; ++j) {
        const Vector preconditioned = apply_right_preconditioner(basis[j]);
        Vector w = jacobian_vector(state, preconditioned, base_residual);
        for (int i = 0; i <= j; ++i) {
            hessenberg(i, j) = w.dot(basis[i]);
            w -= hessenberg(i, j) * basis[i];
        }
        // One reorthogonalization pass keeps the small Krylov basis stable.
        for (int i = 0; i <= j; ++i) {
            const double correction = w.dot(basis[i]);
            hessenberg(i, j) += correction;
            w -= correction * basis[i];
        }
        hessenberg(j + 1, j) = w.norm();
        if (hessenberg(j + 1, j) > 1e-15)
            basis[j + 1] = w / hessenberg(j + 1, j);

        Eigen::VectorXd rhs_small = Eigen::VectorXd::Zero(j + 2);
        rhs_small[0] = beta;
        const Eigen::VectorXd coefficients =
            hessenberg.topLeftCorner(j + 2, j + 1).colPivHouseholderQr().solve(rhs_small);
        solution.setZero();
        for (int i = 0; i <= j; ++i) solution += coefficients[i] * basis[i];
        linear_residual = (rhs_small - hessenberg.topLeftCorner(j + 2, j + 1) * coefficients).norm();
        iterations = j + 1;
        if (!solution.allFinite()) {
            ok = false;
            break;
        }
        if (linear_residual <= 1e-10 * beta)
            break;
        if (hessenberg(j + 1, j) <= 1e-15)
            break;
    }
    ok = ok && solution.allFinite();
    return apply_right_preconditioner(solution);
}

// Newton-Krylov: GMRES on finite-difference Jacobian products, each of which
// pushes all particles, with the mass matrix as right linear preconditioner.
void SimulationEccapfim::solve_fields_newton_krylov(Vector& state, Vector& residual) {
    const double initial_norm = residual.norm();
    bool converged = false;
    for (int outer = 0; outer <= coupled_max_outer_; ++outer) {
        nonlinear_iterations = outer;
        nonlinear_residual = residual.norm();
        if (iteration_log_)
            std::cout << "ECCAPFIM_COUPLED step=" << current_timestep_ << " newton=" << outer
                      << " residual=" << std::setprecision(17) << nonlinear_residual << std::endl;
        if (!residual.allFinite() || !std::isfinite(nonlinear_residual))
            throw std::runtime_error("ECCAPFIM: nonfinite coupled residual");
        if (nonlinear_residual <= coupled_atol_ + coupled_rtol_ * initial_norm) {
            converged = true;
            break;
        }
        if (outer == coupled_max_outer_)
            break;

        if (mass_matrix_every_iteration_)
            mass_matrix_ready_ = false;
        int linear_iterations = 0;
        bool linear_ok = false;
        const Vector direction = gmres(state, -residual, linear_iterations, linear_ok);
        if (!linear_ok || !direction.allFinite())
            throw std::runtime_error("ECCAPFIM: coupled GMRES failed");
        if (iteration_log_)
            std::cout << "ECCAPFIM_COUPLED_GMRES step=" << current_timestep_ << " newton=" << outer
                      << " iterations=" << linear_iterations << std::endl;

        double alpha = 1.0;
        Vector trial_state, trial_residual;
        bool accepted = false;
        for (int line_search = 0; line_search < 10; ++line_search) {
            trial_state = state + alpha * direction;
            trial_residual = coupled_residual(trial_state);
            if (trial_residual.norm() < (1.0 - 1e-4 * alpha) * residual.norm()) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted)
            throw std::runtime_error("ECCAPFIM: coupled Newton line search failed");
        state = trial_state;
        residual = trial_residual;
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: coupled field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
}

// Nonlinear GMRES with right nonlinear preconditioning (Oosterlee & Washio
// 2000; PETSc SNESNGMRES as used by xpic). One preconditioned step
// x^M = x - P^{-1} F(x) uses the frozen-response Jacobian
// P = [[I + dt/2 M, -dt/2 curl_B], [dt/2 curl_E, I]]; then x^A minimizes the
// linearized residual over the history of iterates. Every candidate is judged
// by the full nonlinear residual, i.e. by a complete particle push.
void SimulationEccapfim::solve_fields_ngmres(Vector& state, Vector& residual) {
    constexpr double gamma_a = 2.0, epsilon_b = 0.1, delta_b = 0.9, gamma_c = 2.0;
    constexpr int restart_after = 2;
    const double initial_norm = residual.norm();
    Vector last_evaluated = state;
    std::vector<Vector> history_x, history_f;
    std::vector<double> history_norm;
    int restart_count = 0;
    bool converged = false;
    for (int outer = 0; outer <= coupled_max_outer_; ++outer) {
        nonlinear_iterations = outer;
        nonlinear_residual = residual.norm();
        if (iteration_log_)
            std::cout << "ECCAPFIM_NGMRES step=" << current_timestep_ << " iteration=" << outer
                      << " residual=" << std::setprecision(17) << nonlinear_residual << std::endl;
        if (!residual.allFinite() || !std::isfinite(nonlinear_residual))
            throw std::runtime_error("ECCAPFIM: nonfinite coupled residual");
        if (nonlinear_residual <= coupled_atol_ + coupled_rtol_ * initial_norm) {
            converged = true;
            break;
        }
        if (outer == coupled_max_outer_)
            break;

        history_x.push_back(state);
        history_f.push_back(residual);
        history_norm.push_back(nonlinear_residual);
        if (static_cast<int>(history_x.size()) > coupled_history_) {
            history_x.erase(history_x.begin());
            history_f.erase(history_f.begin());
            history_norm.erase(history_norm.begin());
        }

        // Preconditioned step. fieldBHalf holds B of `state` (last evaluation).
        if (mass_matrix_every_iteration_)
            mass_matrix_ready_ = false;
        if (mass_matrix_check_ && outer == 0)
            check_mass_matrix(state, residual);
        const Vector state_m = state - apply_right_preconditioner(residual);
        const Vector residual_m = coupled_residual(state_m);
        last_evaluated = state_m;
        const double norm_m = residual_m.norm();

        // Least-squares combination over the history, linearized around x^M.
        const int m = static_cast<int>(history_x.size());
        Eigen::MatrixXd dF(residual.size(), m);
        for (int j = 0; j < m; ++j) dF.col(j) = history_f[j] - residual_m;
        const Vector beta = dF.colPivHouseholderQr().solve(-residual_m);
        Vector state_a = state_m;
        for (int j = 0; j < m; ++j) state_a += beta[j] * (history_x[j] - state_m);
        bool select_a = false;
        Vector residual_a;
        double norm_a = std::numeric_limits<double>::infinity();
        const double min_norm = *std::min_element(history_norm.begin(), history_norm.end());
        if (state_a.allFinite() && (state_a - state_m).norm() > 0.0) {
            residual_a = coupled_residual(state_a);
            last_evaluated = state_a;
            norm_a = residual_a.norm();
            if (norm_a * norm_a <= gamma_a * min_norm * min_norm) {
                select_a = true;
            } else {
                double min_distance = std::numeric_limits<double>::infinity();
                for (const auto& x : history_x) min_distance = std::min(min_distance, (state_a - x).norm());
                select_a = epsilon_b * (state_a - state_m).norm() < min_distance ||
                           std::sqrt(norm_a) < delta_b * std::sqrt(min_norm);
            }
            select_a = select_a && std::isfinite(norm_a) && norm_a <= norm_m * 10;
        }
        // Restart the history after repeated growth of the accelerated residual.
        restart_count = (std::isfinite(norm_a) && norm_a > gamma_c * min_norm) ? restart_count + 1 : 0;
        if (restart_count >= restart_after) {
            history_x.clear();
            history_f.clear();
            history_norm.clear();
            restart_count = 0;
        }
        if (select_a) {
            state = state_a;
            residual = residual_a;
        } else {
            state = state_m;
            residual = residual_m;
        }
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: coupled NGMRES field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
    // Particles and sources must belong to the accepted fields.
    if (last_evaluated != state)
        residual = coupled_residual(state);
}

// Diagnostic of the frozen-response approximation (xpic
// -npc_mass_matrix_check): compares P v with the finite-difference Jacobian
// product (F(x + eps v) - F(x)) / eps for random E and B directions.
void SimulationEccapfim::check_mass_matrix(const Vector& state, const Vector& residual) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> u(-1, 1);
    for (int part = 0; part < 2; ++part) {
        Vector v = Vector::Zero(state.size());
        for (int i = 0; i < v.size(); ++i)
            if ((i % 6 < 3) == (part == 0))
                v[i] = u(rng);
        // P v with P = [[maxwell + mass - dt^2/4 curl_B curl_E, -dt/2 curl_B], [dt/2 curl_E, I]].
        ensure_mass_matrix(fieldBHalf + fieldBInit);
        Field3d ve(domain.size(), 3), vb(domain.size(), 3);
        ve.setZero();
        vb.setZero();
        unpack6(v, ve, vb);
        const Eigen::SparseMatrix<double> mass_matrix = mass_only_matrix();
        const Vector pe = pack(ve) + mass_matrix * pack(ve);
        Field3d pe_field(domain.size(), 3);
        pe_field.setZero();
        unpack(pe, pe_field);
        const Field3d p_e = pe_field - (0.5 * dt_) * (maxwell.curlB * vb);
        const Field3d p_b = vb + (0.5 * dt_) * (maxwell.curlE * ve);
        const Vector pv = pack6(p_e, p_b);
        const double eps = 1e-7 * (1.0 + state.norm()) / v.norm();
        const Vector jv = (coupled_residual(state + eps * v) - residual) / eps;
        // The particle part alone: (dt/2) dJ v against the mass matrix.
        Field3d mv_field(domain.size(), 3);
        mv_field.setZero();
        unpack(mass_matrix * pack(ve), mv_field);
        Field3d zero(domain.size(), 3);
        zero.setZero();
        const Vector mv = pack6(mv_field, zero);
        const Vector response = jv - (pv - mv);
        std::cout << "MASS_MATRIX_CHECK step=" << current_timestep_ << " direction=" << (part == 0 ? "E" : "B")
                  << " |Jv|=" << jv.norm() << " |Pv-Jv|/|Jv|=" << (pv - jv).norm() / jv.norm()
                  << " |dt/2 dJ|=" << response.norm()
                  << " |Mv - dt/2 dJ|/|dt/2 dJ|=" << (mv - response).norm() / std::max(response.norm(), 1e-300)
                  << std::endl;
    }
    coupled_residual(state);   // restore particles and sources of `state`
}

// Anderson acceleration of the preconditioned fixed point
// G(x) = x - P^{-1} F(x) (Walker & Ni 2011): one particle push per
// iteration, against two for NGMRES. With history m the update is
// x+ = x + g - (dX + dG) gamma, gamma = argmin |g - dG gamma|, g = G(x) - x.
// Convergence and acceptance use the full nonlinear residual F, and the
// accepted state is always the last evaluated one.
void SimulationEccapfim::solve_fields_anderson(Vector& state, Vector& residual) {
    const double initial_norm = residual.norm();
    std::vector<Vector> history_x, history_g;
    bool converged = false;
    for (int outer = 0; outer <= coupled_max_outer_; ++outer) {
        nonlinear_iterations = outer;
        nonlinear_residual = residual.norm();
        if (iteration_log_)
            std::cout << "ECCAPFIM_ANDERSON step=" << current_timestep_ << " iteration=" << outer
                      << " residual=" << std::setprecision(17) << nonlinear_residual << std::endl;
        if (!residual.allFinite() || !std::isfinite(nonlinear_residual))
            throw std::runtime_error("ECCAPFIM: nonfinite coupled residual");
        if (nonlinear_residual <= coupled_atol_ + coupled_rtol_ * initial_norm) {
            converged = true;
            break;
        }
        if (outer == coupled_max_outer_)
            break;
        if (mass_matrix_every_iteration_)
            mass_matrix_ready_ = false;
        const Vector g = -apply_right_preconditioner(residual);
        Vector next = state + g;
        const int m = static_cast<int>(history_x.size());
        if (m > 0) {
            Eigen::MatrixXd dX(state.size(), m), dG(state.size(), m);
            for (int j = 0; j < m; ++j) {
                dX.col(j) = state - history_x[j];
                dG.col(j) = g - history_g[j];
            }
            const Vector gamma = dG.colPivHouseholderQr().solve(g);
            const Vector accelerated = next - (dX + dG) * gamma;
            if (accelerated.allFinite())
                next = accelerated;
        }
        history_x.push_back(state);
        history_g.push_back(g);
        if (static_cast<int>(history_x.size()) > coupled_history_) {
            history_x.erase(history_x.begin());
            history_g.erase(history_g.begin());
        }
        state = next;
        residual = coupled_residual(state);
    }
    if (!converged)
        throw std::runtime_error("ECCAPFIM: coupled Anderson field solve failed, residual=" +
                                 std::to_string(nonlinear_residual));
}

void SimulationEccapfim::solve_fields_coupled() {
    RECORD_TIMER;
    Vector state = pack6(fieldE, fieldB);
    Vector residual = coupled_residual(state);
    if (coupled_solver_ == CoupledSolver::Anderson)
        solve_fields_anderson(state, residual);
    else if (coupled_solver_ == CoupledSolver::Ngmres)
        solve_fields_ngmres(state, residual);
    else
        solve_fields_newton_krylov(state, residual);
    // B_h is the midpoint field; promote both fields to integer time levels.
    unpack6(state, fieldEp, fieldBHalf);
    fieldEn = 2.0 * fieldEp - fieldE;
    fieldBn = 2.0 * fieldBHalf - fieldB;
    algorithmsECCAPFIM::copy_periodic_nodes(fieldEn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    algorithmsECCAPFIM::copy_periodic_nodes(fieldBn, domain, &bc_handler);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
}
