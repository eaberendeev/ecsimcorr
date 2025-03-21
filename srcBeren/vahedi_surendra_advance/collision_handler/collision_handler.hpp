// collision_handler.hpp

#ifndef COLLISION_HANDLER_HPP
#define COLLISION_HANDLER_HPP

//#include "clsVec.hpp"
#include "cross_section.hpp"


// Перечисление типов столкновений
enum class CollisionType {
    IONIZATION,
    CHARGE_EXCHANGE,
    NULL_COLLISION
};

// Функция вычисления вероятности столкновения
double compute_collision_probability(double freq);

// Функция проверки, произошло ли столкновение
bool check_collision(double P_collision);

// Функция выбора типа столкновения
CollisionType select_collision_type(double E, double mcp, double nn, double freq_max);

#endif // COLLISION_HANDLER_HPP
