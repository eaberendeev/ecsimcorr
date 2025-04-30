// collision_handler.cpp

// #include "collision_handler.hpp"

#include <random>

#include "collision_utils.h"
#include "collisions_with_neutrals.h"
#include "cross_section.hpp"

using namespace std;


// Функция проверки, произошло ли столкновение
bool ColliderWithNeutrals::check_collision(double P_collision) {
    // Генератор случайных чисел для вероятностей
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);


    return dist(gen) < P_collision;
}

// Функция выбора типа столкновения
CollisionType ColliderWithNeutrals::select_collision_type(double E, double mcp,
                                                          double nn,
                                                          double freq_max) {
    double v_mod = compute_velocity(E, mcp);

    // Генератор случайных чисел для вероятностей
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    
    double r = dist(gen) * freq_max;
    if (mcp == 1.) {
        double freq_e = v_mod * nn * Sigma_e(E);
       // std::cout << "freq_e: " << freq_e << std::endl;
        return (r < freq_e) ? CollisionType::IONIZATION : CollisionType::NULL_COLLISION;
    } else {
        double freq_p = v_mod * nn * Sigma_p(E);
        double freq_cx = v_mod * nn * Sigma_cx(E);
        //std::cout << "freq_p: " << freq_p << std::endl;
        //std::cout << "freq_cx: " << freq_cx << std::endl;
        if (r < freq_p) {
            return CollisionType::NULL_COLLISION; // IONIZATION;   // NULL_COLLISION
        } else if (r < freq_p + freq_cx) {
            return CollisionType::CHARGE_EXCHANGE;
        } else {
            return CollisionType::NULL_COLLISION;
        }
    }
}

// for (celL : cells){
//     v_max = max(v_cell);
//     max_sigma = max(sigma_cell);

// }