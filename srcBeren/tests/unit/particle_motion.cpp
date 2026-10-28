// Single particle in prescribed constant fields: no field solver, sources or walls.
// Velocities are in c, time in 1/omega_pe, E/B in the code normalization.
#include <iomanip>

#include "boris_pusher.h"
#include "test_support.h"

namespace {
struct PhasePoint {
    Vector3R position;
    Vector3R velocity;
};

// Continuous Lorentz equation: dv/dt = (q/m)*(E + v cross B).
PhasePoint exact_orbit(const PhasePoint& initial, double qm, const Vector3R& electric, const Vector3R& magnetic,
                       double time) {
    if (magnetic.norm() == 0) {
        const Vector3R acceleration = qm * electric;
        return {initial.position + time * initial.velocity + (0.5 * time * time) * acceleration,
                initial.velocity + time * acceleration};
    }
    const Vector3R axis = magnetic / magnetic.norm();
    const double omega = qm * magnetic.norm();
    const Vector3R drift = electric.cross(magnetic) / magnetic.squaredNorm();
    const Vector3R parallel_velocity = initial.velocity.dot(axis) * axis;
    const Vector3R parallel_acceleration = (qm * electric.dot(axis)) * axis;
    const Vector3R rotating_velocity = initial.velocity - parallel_velocity - drift;
    const Vector3R quadrature = rotating_velocity.cross(axis);
    const double angle = omega * time;
    return {initial.position + time * (drift + parallel_velocity) + (0.5 * time * time) * parallel_acceleration +
                (std::sin(angle) / omega) * rotating_velocity + ((1 - std::cos(angle)) / omega) * quadrature,
            drift + parallel_velocity + time * parallel_acceleration + std::cos(angle) * rotating_velocity +
                std::sin(angle) * quadrature};
}

struct OrbitError {
    double position;
    double velocity;
};

OrbitError track(const std::string& name, bool staggered, double dt, const Vector3R& electric,
                 const Vector3R& magnetic) {
    constexpr double qm = -0.7;
    constexpr double final_time = 12;
    const PhasePoint initial{Vector3R(1, -2, 0.5), Vector3R(0.2, -0.1, 0.05)};
    const PhasePoint coordinate_level = staggered ? exact_orbit(initial, qm, electric, magnetic, -0.5 * dt) : initial;
    Particle particle(coordinate_level.position, initial.velocity);
    const int steps = std::lround(final_time / dt);
    for (int step = 0; step < steps; ++step) {
        // ECSIM: x_(n-1/2) -> x_(n+1/2), then v_n -> v_(n+1).
        // Corr in uniform fields: two half drifts around a full Boris update.
        particle.move(staggered ? dt : 0.5 * dt);
        borisPusher::update_vEB(particle, qm, electric, magnetic, dt);
        if (!staggered)
            particle.move(0.5 * dt);
        if (electric.norm() == 0)
            check::near(particle.velocity.squaredNorm(), initial.velocity.squaredNorm(), 1e-13,
                        "constant B must conserve kinetic energy at every step");
    }
    const auto exact_velocity = exact_orbit(initial, qm, electric, magnetic, final_time);
    const auto exact_position = exact_orbit(initial, qm, electric, magnetic, final_time - (staggered ? 0.5 * dt : 0));
    const OrbitError error{(particle.coord - exact_position.position).norm(),
                           (particle.velocity - exact_velocity.velocity).norm()};
    check::require(std::isfinite(error.position) && std::isfinite(error.velocity), "nonfinite particle orbit");
    std::cout << name << (staggered ? " ECSIM" : " Corr") << " dt=" << dt << " |x-x_exact|=" << error.position
              << " |v-v_exact|=" << error.velocity << '\n';
    return error;
}

void constant_fields() {
    for (bool staggered : {false, true}) {
        const auto electric = track("constant E, B=0", staggered, 0.2, Vector3R(0.01, -0.02, 0.03), Vector3R(0.0));
        check::near(electric.position, 0, 1e-12, "constant electric acceleration: analytic position");
        check::near(electric.velocity, 0, 1e-13, "constant electric acceleration: analytic velocity");
        for (bool crossed : {false, true}) {
            const std::string name =
                crossed ? "constant E+B (gyro + drift + parallel acceleration)" : "constant B (gyro)";
            const Vector3R e = crossed ? Vector3R(0.01, 0.02, -0.005) : Vector3R(0.0);
            const Vector3R b(0.2, -0.3, 0.5);
            const auto coarse = track(name, staggered, 0.2, e, b);
            const auto fine = track(name, staggered, 0.1, e, b);
            // Boris and the centered drifts have second order global error:
            // halving dt should reduce the orbit error by approximately four.
            check::require(coarse.position > 1e-8 && coarse.velocity > 1e-8,
                           "orbit must resolve a nonzero phase error");
            const double x_ratio = coarse.position / fine.position, v_ratio = coarse.velocity / fine.velocity;
            std::cout << "convergence ratios: position=" << x_ratio << " velocity=" << v_ratio
                      << " expected [3.5,4.5]\n";
            check::require(std::isfinite(x_ratio) && x_ratio >= 3.5 && x_ratio <= 4.5, "position is not second order");
            check::require(std::isfinite(v_ratio) && v_ratio >= 3.5 && v_ratio <= 4.5, "velocity is not second order");
            check::near(fine.position, 0, 0.002, "absolute spatial orbit accuracy at dt=0.1");
            check::near(fine.velocity, 0, 0.001, "absolute velocity orbit accuracy at dt=0.1");
        }
    }
}

void pusher() {
    const Vector3R initial(0.2, -0.1, 0.4), electric(0.03, -0.02, 0.01), magnetic(0, 0, 0.7);
    Particle p(Vector3R(0, 0, 0), initial);
    borisPusher::update_vEB(p, -0.5, electric, Vector3R(0, 0, 0), 0.4);
    check::vector(p.velocity, initial - 0.2 * electric, 1e-15, "Boris electric impulse");
    p.velocity = Vector3R(1, 0, 0);
    borisPusher::update_vEB(p, 1, Vector3R(0, 0, 0), magnetic, 0.1);
    check::require(p.velocity.y() < 0, "Lorentz rotation has wrong sign");
    p.velocity = initial;
    for (int step = 0; step < 2000; ++step) borisPusher::update_vEB(p, -0.5, Vector3R(0, 0, 0), magnetic, 0.1);
    check::near(p.velocity.squaredNorm(), initial.squaredNorm(), 1e-13, "magnetic rotation changed kinetic energy");
    p.velocity = initial;
    borisPusher::update_vEB(p, -0.5, electric, magnetic, 0.4);
    borisPusher::update_vEB(p, -0.5, electric, magnetic, -0.4);
    check::vector(p.velocity, initial, 1e-14, "Boris reversibility");
}

}   // namespace

int main() {
    return check::run([] {
        std::cout << std::scientific << std::setprecision(8);
        constant_fields();
        pusher();   // Electric impulse, rotation sign and time reversibility.
        std::cout << "PASS: particle_motion (analytic trajectories, order, energy, reversibility)\n";
    });
}
