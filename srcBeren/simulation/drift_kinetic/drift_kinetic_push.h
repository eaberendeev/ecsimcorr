#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include "core/vector3.h"

namespace driftKinetic {

// Guiding-center state. `mu` is the magnetic moment (energy per unit B); the
// parallel velocity and position are the implicit midpoint unknowns.
struct GuidingCenter {
    Vector3R position;
    double v_parallel = 0.0;
    double mu = 0.0;
    double charge = 0.0;
    double mass = 1.0;
    double weight = 1.0;
};

// Fields interpolated at the orbit midpoint by the caller.
struct Fields {
    Vector3R electric;
    double lenB = 0.0;
    Vector3R b;   // unit magnetic direction
    Vector3R gradB;
    Vector3R rotB;
};

using FieldFunction = std::function<Fields(const Vector3R& start, const Vector3R& end)>;

struct PushResult {
    GuidingCenter state;
    int iterations = 0;
    bool converged = false;
    double position_residual = 0.0;
    double velocity_residual = 0.0;
};

// Fully implicit guiding-center pusher, ported from xpic DriftKineticPush.
// Unknowns: r^{n+1} and v_par^{n+1}; the equations use V_par = (v_par^n + v_par^{n+1}) / 2
// and the fields of the orbit r^n -> r^{n+1} returned by `fields`.
// The sequential Picard/Aitken map updates r and v_parallel together; the
// original midpoint equations are used for acceptance and convergence.
class GuidingCenterPush {
   public:
    void set_tolerances(double atol, double rtol, int max_iterations) {
        if (!(atol > 0.0) || !(rtol >= 0.0) || max_iterations < 1 || !std::isfinite(atol) || !std::isfinite(rtol))
            throw std::invalid_argument("drift-kinetic tolerances must be finite and positive");
        atol_ = atol;
        rtol_ = rtol;
        max_iterations_ = max_iterations;
    }

    PushResult push(const GuidingCenter& initial, double dt, const FieldFunction& fields,
                    bool use_aitken = true) const {
        return push(initial, initial, dt, fields, use_aitken);
    }

    // `initial` is the state at the start of the (sub)step, `guess` the
    // starting iterate for the end state; the outer field iteration passes the
    // orbit of its previous evaluation (xpic warm start).
    PushResult push(const GuidingCenter& initial, const GuidingCenter& guess, double dt, const FieldFunction& fields,
                    bool use_aitken) const {
        if (!(dt > 0.0) || !std::isfinite(dt))
            throw std::invalid_argument("drift-kinetic push requires a finite positive dt");
        if (!(initial.mass > 0.0) || !std::isfinite(initial.mass) || initial.charge == 0.0 ||
            !std::isfinite(initial.charge))
            throw std::invalid_argument("drift-kinetic push requires finite nonzero charge and positive mass");
        if (!fields)
            throw std::invalid_argument("drift-kinetic push requires a field callback");

        const double qm = initial.charge / initial.mass;
        const double mp = initial.mass;
        // `mu` is energy per unit B, so the mirror force is -mu grad|B|.
        const double mu_p = initial.mu;
        GuidingCenter pn = initial;
        pn.position = guess.position;
        pn.v_parallel = guess.v_parallel;
        const GuidingCenter p0 = initial;

        Fields f;
        Vector3R B_eff, b_eff, F_eff, V_dgc;
        double lenB_eff = 0.0, V_par = initial.v_parallel, a_par = 0.0;
        Vector3R V_dgc_prev;
        double a_par_prev = 0.0;
        bool has_prev = false;

        const auto finite = [](const Vector3R& value) {
            return std::isfinite(value.x()) && std::isfinite(value.y()) && std::isfinite(value.z());
        };

        const auto update_drift = [&]() {
            B_eff = f.lenB * f.b + (V_par / qm) * f.rotB;
            lenB_eff = f.b.dot(B_eff);
            if (!(std::abs(lenB_eff) > 0.0) || !std::isfinite(lenB_eff))
                throw std::runtime_error("drift-kinetic effective magnetic field is degenerate");
            b_eff = B_eff / lenB_eff;
            V_dgc = V_par * b_eff + F_eff.cross(f.b) / (initial.charge * lenB_eff);
        };

        const auto evaluate = [&](const GuidingCenter& state) {
            V_par = 0.5 * (state.v_parallel + p0.v_parallel);
            f = fields(p0.position, state.position);
            if (!(f.lenB > 0.0) || !std::isfinite(f.lenB))
                throw std::runtime_error("drift-kinetic requires a finite positive magnetic field");
            if (!finite(f.b) || !finite(f.gradB) || !finite(f.rotB) || !finite(f.electric))
                throw std::runtime_error("drift-kinetic field callback returned non-finite data");
            f.b = f.b.normalized();
            f.rotB = (f.b.cross(f.gradB) + f.rotB) / f.lenB;
            F_eff = initial.charge * f.electric - mu_p * f.gradB;
            update_drift();
            a_par = F_eff.dot(b_eff) / mp;
            if (!finite(V_dgc) || !std::isfinite(a_par))
                throw std::runtime_error("drift-kinetic RHS is non-finite");
        };

        auto position_residual = [&]() { return (pn.position - p0.position - dt * V_dgc).length(); };
        auto velocity_residual = [&]() { return std::abs(pn.v_parallel - p0.v_parallel - dt * a_par); };
        auto residual_ok = [&](double FRk, double FVk) {
            const double r_scale = dt * atol_ + rtol_ * (dt * V_dgc.length());
            const double v_scale = atol_ + rtol_ * (std::abs(V_par) + std::abs(dt * a_par));
            return FRk <= r_scale && FVk <= v_scale;
        };
        auto absolute_ok = [&]() { return residual_ok(position_residual(), velocity_residual()); };
        auto previous_ok = [&]() {
            return residual_ok(dt * (V_dgc_prev - V_dgc).length(), dt * std::abs(a_par_prev - a_par));
        };
        auto increment_ok = [&](const std::array<double, 4>& increment) {
            const Vector3R dr(increment[0], increment[1], increment[2]);
            return residual_ok(dr.length(), std::abs(increment[3]) / dt);
        };

        PushResult result;
        result.state = pn;
        int iterations = 0;
        bool converged = false;
        double FRk = 0.0, FVk = 0.0;

        if (!use_aitken) {
            for (iterations = 0; iterations < max_iterations_; ++iterations) {
                evaluate(pn);
                const bool ok = has_prev ? previous_ok() : absolute_ok();
                if (ok) {
                    converged = true;
                    break;
                }
                FRk = has_prev ? dt * (V_dgc_prev - V_dgc).length() : position_residual();
                FVk = has_prev ? dt * std::abs(a_par_prev - a_par) : velocity_residual();
                pn.position = p0.position + dt * V_dgc;
                pn.v_parallel = p0.v_parallel + dt * a_par;
                V_dgc_prev = V_dgc;
                a_par_prev = a_par;
                has_prev = true;
            }
            if (!converged) {
                evaluate(pn);
                FRk = position_residual();
                FVk = velocity_residual();
                converged = residual_ok(FRk, FVk);
            }
        } else {
            Vector3R dr = pn.position - p0.position;
            double dv = pn.v_parallel - p0.v_parallel;
            std::array<double, 4> previous{};
            double omega = 1.0;
            bool has_history = false;
            iterations = 0;
            evaluate(pn);

            for (;;) {
                const Vector3R rr = dt * V_dgc - dr;
                FRk = rr.length();
                FVk = std::abs(dt * (dt * a_par - dv)) / dt;
                converged = increment_ok({rr.x(), rr.y(), rr.z(), dt * (dt * a_par - dv)});
                if (converged || iterations >= 4 * max_iterations_)
                    break;

                double norm2 = FRk * FRk + std::pow(FVk * dt, 2);
                if (!std::isfinite(norm2))
                    break;

                const double dv_next = dt * a_par;
                V_par = 0.5 * (p0.v_parallel + (p0.v_parallel + dv_next));
                update_drift();
                const Vector3R next = dt * V_dgc - dr;
                const std::array<double, 4> f_increment{next.x(), next.y(), next.z(), dt * (dv_next - dv)};

                if (has_history) {
                    double numerator = 0.0, denominator = 0.0, scale = 0.0;
                    for (int i = 0; i < 4; ++i) {
                        const double difference = f_increment[i] - previous[i];
                        numerator += previous[i] * difference;
                        denominator += difference * difference;
                        scale += f_increment[i] * f_increment[i] + previous[i] * previous[i];
                    }
                    if (denominator > std::numeric_limits<double>::epsilon() * scale) {
                        const double candidate = -omega * numerator / denominator;
                        if (std::isfinite(candidate))
                            omega = std::clamp(candidate, 0.05, 1.0);
                    }
                }

                const GuidingCenter current = pn;
                for (;;) {
                    const Vector3R trial_dr = dr + omega * Vector3R(f_increment[0], f_increment[1], f_increment[2]);
                    const double trial_dv = dv + omega * f_increment[3] / dt;
                    pn.position = p0.position + trial_dr;
                    pn.v_parallel = p0.v_parallel + trial_dv;
                    ++iterations;

                    double trial_norm2 = std::numeric_limits<double>::max();
                    if (finite(pn.position) && std::isfinite(pn.v_parallel)) {
                        evaluate(pn);
                        const Vector3R trial_rr = dt * V_dgc - trial_dr;
                        const double trial_rv = dt * (dt * a_par - trial_dv);
                        FRk = trial_rr.length();
                        FVk = std::abs(trial_rv) / dt;
                        trial_norm2 = trial_rr.squaredNorm() + trial_rv * trial_rv;
                        converged = increment_ok({trial_rr.x(), trial_rr.y(), trial_rr.z(), trial_rv});
                    }

                    if (std::isfinite(trial_norm2) && (converged || trial_norm2 <= 4.0 * norm2)) {
                        dr = trial_dr;
                        dv = trial_dv;
                        previous = f_increment;
                        has_history = true;
                        break;
                    }
                    if (omega <= 0.05 || iterations >= 4 * max_iterations_) {
                        pn = current;
                        evaluate(pn);
                        const Vector3R restore_rr = dt * V_dgc - dr;
                        FRk = restore_rr.length();
                        FVk = std::abs(dt * a_par - dv);
                        converged =
                            increment_ok({restore_rr.x(), restore_rr.y(), restore_rr.z(), dt * (dt * a_par - dv)});
                        result.state = pn;
                        result.iterations = iterations;
                        result.converged = converged;
                        result.position_residual = FRk;
                        result.velocity_residual = FVk;
                        return result;
                    }
                    omega = std::max(0.05, 0.5 * omega);
                }
            }
        }

        result.state = pn;
        result.iterations = iterations;
        result.converged = converged;
        result.position_residual = FRk;
        result.velocity_residual = FVk;
        return result;
    }

   private:
    double atol_ = 1e-16;
    double rtol_ = 1e-13;
    int max_iterations_ = 100;
};

}   // namespace driftKinetic
