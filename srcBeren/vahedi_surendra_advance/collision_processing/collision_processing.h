/**
 * @file collision_processing.h
 * @author Morozov O. P.
 * @brief Declares collision processing helpers for particle-neutral interactions.
 */
#ifndef COLLISION_PROCESSING_HPP
#define COLLISION_PROCESSING_HPP

#include <tuple>

#include "Vec.h"
#include "cross_section.h"

double3 get_scattered_velocity(double speed);

double3 get_electron_scattered_velocity(double3 velocity, double energy);

double3 get_proton_scattered_velocity(double3 velocity);

std::tuple<bool, double3, double3> process_collision(
    CollisionType collision_type, bool is_electron,
    double3& vcp, double3& vn,
    double mcp, double mn
);

#endif // COLLISION_PROCESSING_HPP
