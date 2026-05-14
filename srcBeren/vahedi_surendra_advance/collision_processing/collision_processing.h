#ifndef COLLISION_PROCESSING_HPP
#define COLLISION_PROCESSING_HPP

#include <tuple>

#include "cross_section.h"
#include "vector3.h"

Vector3R get_scattered_velocity(double speed);

Vector3R get_electron_scattered_velocity(Vector3R velocity, double energy);

Vector3R get_proton_scattered_velocity(Vector3R velocity);

std::tuple<bool, Vector3R, Vector3R> process_collision(CollisionType collision_type, bool is_electron, Vector3R& vcp,
                                                       Vector3R& vn, double mcp, double mn);

#endif   // COLLISION_PROCESSING_HPP
