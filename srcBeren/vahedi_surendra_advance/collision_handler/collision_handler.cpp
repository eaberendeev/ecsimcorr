// collision_handler.cpp

#include "collision_handler.hpp"
#include "cross_section.hpp"
#include "utils.hpp"

#include <random>

using namespace std;

// Функция вычисления вероятности столкновения
// P = 1 - exp(-sigma * n * v * dt)
double compute_collision_probability(double freq) {
    return 1.0 - exp(- freq * dt);
}

// Функция проверки, произошло ли столкновение
bool check_collision(double P_collision) {
    // Генератор случайных чисел для вероятностей
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);


    return dist(gen) < P_collision;
}

// Функция выбора типа столкновения
CollisionType select_collision_type(double E, double mcp, double nn, double freq_max) {
    double v_mod = compute_velocity(E, mcp);

    // Генератор случайных чисел для вероятностей
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    
    double r = dist(gen) * freq_max;
    if (mcp == 1.) {
        double freq_e = v_mod * nn * Sigma_e(E);
        return (r < freq_e) ? CollisionType::IONIZATION : CollisionType::NULL_COLLISION;
    } else {
        double freq_p = v_mod * nn * Sigma_p(E);
        double freq_cx = v_mod * nn * Sigma_cx(E);
        if (r < freq_p) {
            return CollisionType::IONIZATION; // 
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