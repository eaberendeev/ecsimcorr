// collision_processing.cpp

#include "collision_processing.hpp"
#include "utils.hpp"
#include <cmath>
#include <random>
#include "Vec.h"

using namespace std;

// Генератор случайных чисел
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> dist(0.0, 1.0);


// Функция генерации случайного рассеяния (Vahedi)
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

// Функция генерации угла рассеяния для электрона (формула 9 из статьи Vahedi)
double3 get_electron_scattered_velocity(double3 velocity, double energy) {
    double g = energy; // Нормированная энергия
    double cos_theta = (2. + g - 2. * pow(1. + g, dist(gen))) / g;
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double phi = 2 * M_PI * dist(gen);
    
    double speed = velocity.norm();
    return double3(
        speed * sin_theta * cos(phi),
        speed * sin_theta * sin(phi),
        speed * cos_theta
    );
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


// Обработка столкновения и обновление скоростей частиц
tuple<bool, double3, double3> process_collision(
    CollisionType collision_type,
    double3& vcp, double3& vn,
    double mcp, double mn
) {
    switch (collision_type) {
        case CollisionType::IONIZATION: {
            if (mcp == 1.0) { // Электронная ионизация (модель Опала)
                double R = dist(gen);
                double E_inc = compute_energy(vcp - vn, mcp);
                double E_ej = B_param * tan(R * atan((E_inc - E_ion) / (2.0 * B_param)));
                double E_scat = E_inc - E_ion - E_ej;

                cout << squaredNorm(vcp - vn) << endl;

                double3 vcp_scat = get_electron_scattered_velocity((vcp - vn), E_inc);
                cout << squaredNorm(vcp_scat) << endl;
                vcp = vcp_scat + vn;

                double3 new_electron_velocity = vn + get_scattered_velocity(compute_velocity(E_ej, mcp));
                //cout << "vcp2: " << norm(vcp) + norm(ejected_electron_velocity) << endl; 
                
                double3 new_ion_velocity = vn;
                
                return {true, new_electron_velocity, new_ion_velocity};
            } else { // Протонная ионизация (модель Vahedi)
                double E_inc = compute_energy(vcp - vn, mcp);

                vcp = (vcp - vn) * sqrt(1. - E_ion/E_inc) + vn;

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
            swap(vcp, vn);
            return {false, vcp, vn};
        }
        case CollisionType::NULL_COLLISION:
        default:
            return {false, vcp, vn};
    }
}
