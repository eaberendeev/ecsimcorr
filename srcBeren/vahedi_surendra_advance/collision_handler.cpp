/**
 * @file collision_handler.cpp
 * @author Morozov O. P.
 * @brief Implements stochastic collision checks and type selection logic.
 */
#include <random>

#include "collision_utils.h"
#include "collisions_with_neutrals.h"

using namespace std;


// Функция проверки, произошло ли столкновение
bool ColliderWithNeutrals::check_collision(double P_collision) {
    if (P_collision <= 1.e-16)
        return false;
    // если P >= 1.0, сразу true (без генерации)
    if (P_collision >= 1.0)
        return true;
    double u = rng_uniform01();
    return u < P_collision;
}

// Функция выбора типа столкновения
CollisionType ColliderWithNeutrals::select_collision_type(bool is_electron, double ion_freq,
                                                          double cx_freq,
                                                          double freq_bound) {
    if (freq_bound <= 0.0) {
        return CollisionType::NULL_COLLISION;
    }

    // Генератор случайных чисел для вероятностей
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    double r = dist(gen) * freq_bound;
    if (r < ion_freq) {
        return CollisionType::IONIZATION;
    }
    if (r < ion_freq + cx_freq && !is_electron) {
        return CollisionType::CHARGE_EXCHANGE;
    }
    return CollisionType::NULL_COLLISION;
}
