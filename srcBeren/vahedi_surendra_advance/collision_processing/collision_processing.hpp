// collision_processing.hpp

#ifndef COLLISION_PROCESSING_HPP
#define COLLISION_PROCESSING_HPP

#include "Vec.h"
//#include "collision_handler.hpp"
#include "cross_section.hpp"
#include <tuple>

double3 get_scattered_velocity(double speed);


// Функция генерации угла рассеяния для электрона (формула 9 из статьи Vahedi)
double3 get_electron_scattered_velocity(double3 velocity, double energy);

// Функция обработки столкновения
std::tuple<bool, double3, double3> process_collision(
    CollisionType collision_type,
    double3& vcp, double3& vn,
    double mcp, double mn
);

#endif // COLLISION_PROCESSING_HPP
