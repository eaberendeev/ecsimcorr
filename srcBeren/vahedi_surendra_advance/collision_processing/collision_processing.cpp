/**
 * @file collision_processing.cpp
 * @author Morozov O. P.
 * @brief Implements routines that apply specific collision outcomes to particles.
 */
#include "collision_processing.h"
#include "collision_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

using namespace std;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> dist(0.0, 1.0);



double3 get_scattered_velocity(double speed) {
    double cos_theta = 1 - 2 * dist(gen);
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double phi = 2 * M_PI * dist(gen);
    
    return double3(
        speed * sin_theta * cos(phi),
        speed * sin_theta * sin(phi),
        speed * cos_theta
    );
}

double3 get_electron_scattered_velocity(double3 velocity, double energy) {
    const double g = std::max(energy, std::numeric_limits<double>::epsilon());
    double cos_theta = (2. + g - 2. * std::pow(1. + g, dist(gen))) / g;
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);
    const double sin_theta = std::sqrt(std::max(0.0, 1. - cos_theta * cos_theta));
    const double phi = 2. * M_PI * dist(gen);

    double initial_norm = velocity.norm();
    if (initial_norm == 0.0) {
        return get_scattered_velocity(1.0);
    }

    double3 direction = velocity / initial_norm;
    double3 reference(1.0, 0.0, 0.0);
    double3 perpendicular = cross(direction, reference);
    double perpendicular_norm = perpendicular.norm();
    if (perpendicular_norm < 1e-12) {
        reference = double3(0.0, 1.0, 0.0);
        perpendicular = cross(direction, reference);
        perpendicular_norm = perpendicular.norm();
        if (perpendicular_norm < 1e-12) {
            reference = double3(0.0, 0.0, 1.0);
            perpendicular = cross(direction, reference);
            perpendicular_norm = perpendicular.norm();
            if (perpendicular_norm < 1e-12) {
                return get_scattered_velocity(1.0);
            }
        }
    }

    perpendicular = perpendicular / perpendicular_norm;
    double3 secondary = cross(direction, perpendicular);
    const double secondary_norm = secondary.norm();
    if (secondary_norm > 0.0) {
        secondary = secondary / secondary_norm;
    }

    return sin_theta * (std::cos(phi) * perpendicular + std::sin(phi) * secondary) +
           cos_theta * direction;
}

double3 get_proton_scattered_velocity(double3 velocity) {
    double cos_hi = sqrt(1 - dist(gen));
    double sin_hi = sqrt(1 - cos_hi * cos_hi);
    double phi = 2 * M_PI * dist(gen);

    double speed = velocity.norm();
    
    return double3(
        speed * sin_hi * cos(phi),
        speed * sin_hi * sin(phi),
        speed * cos_hi
    );
}


namespace {

double sample_secondary_energy_opal(double available_energy) {
    if (available_energy <= 0.0) {
        return 0.0;
    }
    const double angle_limit = std::atan(available_energy / (2.0 * B_param));
    if (angle_limit <= 0.0) {
        return 0.0;
    }
    const double draw = dist(gen);
    const double candidate = B_param * std::tan(draw * angle_limit);
    if (!std::isfinite(candidate)) {
        return available_energy;
    }
    return std::clamp(candidate, 0.0, available_energy);
}

double sample_secondary_energy_approx(double available_energy) {
    if (available_energy <= 0.0) {
        return 0.0;
    }
    return dist(gen) * available_energy;
}

double3 normalize_direction(double3 vec) {
    const double magnitude = vec.norm();
    if (magnitude <= 0.0) {
        return get_scattered_velocity(1.0);
    }
    return vec / magnitude;
}

}

// Обработка столкновения и обновление скоростей частиц
tuple<bool, double3, double3> process_collision(
    CollisionType collision_type, bool is_electron,
    double3& vcp, double3& vn,
    double mcp, double mn
) {
    switch (collision_type) {
        case CollisionType::IONIZATION: {
            if (is_electron) {
                const double E_inc = compute_energy(vcp - vn, mcp);
                if (E_inc <= E_ion) {
                    return {false, vcp, vn};
                }

                const double available_energy = E_inc - E_ion;
                const auto model = get_electron_ionization_model();
                double E_ej = (model == ElectronIonizationModel::Opal)
                                  ? sample_secondary_energy_opal(available_energy)
                                  : sample_secondary_energy_approx(available_energy);
                E_ej = std::clamp(E_ej, 0.0, available_energy);
                const double E_scat = available_energy - E_ej;

                const double scatter_speed = compute_velocity(E_scat, mcp);
                double3 scatter_direction = normalize_direction(
                    get_electron_scattered_velocity((vcp - vn), E_inc));
                vcp = scatter_direction * scatter_speed + vn;

                const double ejected_speed = compute_velocity(E_ej, mcp);
                double3 new_electron_velocity = vn + get_scattered_velocity(ejected_speed);
                double3 new_ion_velocity = vn;

                return {true, new_electron_velocity, new_ion_velocity};
            } else {
                double E_inc = compute_energy(vcp - vn, mcp);
                if (E_inc <= E_ion) {
                    return {false, vcp, vn};
                }

                const double energy_ratio = std::max(0.0, 1.0 - E_ion / E_inc);
                vcp = (vcp - vn) * std::sqrt(energy_ratio) + vn;

                double3 v_com = v_center_of_mass(vcp, vn, mcp, mn);
                double3 vcp_com = vcp - v_com;
                double3 v_scat = get_proton_scattered_velocity(vcp - v_com);

                vcp = v_scat + v_com;
                double3 new_ion_velocity = (vn + mcp/mn * (vcp_com - v_scat));
                double3 new_electron_velocity = new_ion_velocity;
                
                return {true, new_electron_velocity, new_ion_velocity};
            }
        }
        case CollisionType::CHARGE_EXCHANGE: {
            if(!is_electron){
                swap(vcp, vn);
                return {false, vcp, vn};
            } else {
                return {false, vcp, vn};
            }
        }
        case CollisionType::NULL_COLLISION:
        default:
            return {false, vcp, vn};
    }
}
