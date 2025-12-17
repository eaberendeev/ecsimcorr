#ifndef VOXEL_TRAVERSAL_H
#define VOXEL_TRAVERSAL_H

#include <vector>
#include "vector3.h"

/**
 * @brief returns all the voxels that are traversed by a ray going from start to
 * end
 * @param start : continous world position where the ray starts
 * @param end   : continous world position where the ray end
 * @return vector of voxel ids hit by the ray in temporal order
 *
 * theory: J. Amanatides, A. Woo. A Fast Voxel Traversal Algorithm for Ray
 * Tracing. Eurographics '87
 *
 * realization from: https://github.com/francisengelmann/fast_voxel_traversal
 */
std::vector<Vector3I> voxel_traversal(const Vector3R& ray_start,
                                  const Vector3R& ray_end,
                                  const double bin_size);

double find_ray_voxel_intersection_parameter(const Vector3R& ray_start,
                                             const Vector3R& ray_end,
                                             const Vector3I& current_voxel,
                                             const Vector3I& next_voxel,
                                             double bin_size);
Vector3R get_point_in_ray(const Vector3R& ray_start, const Vector3R& ray_end,
                         const double t);
                         
#endif // VOXEL_TRAVERSAL_H
