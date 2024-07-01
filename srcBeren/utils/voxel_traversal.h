#ifndef VOXEL_TRAVERSAL_H
#define VOXEL_TRAVERSAL_H

#include <vector>
#include "Vec.h"

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
std::vector<int3> voxel_traversal(const double3& ray_start,
                                  const double3& ray_end,
                                  const double bin_size);

double find_ray_voxel_intersection_parameter(const double3& ray_start,
                                             const double3& ray_end,
                                             const int3& current_voxel,
                                             const int3& next_voxel,
                                             double bin_size);
double3 get_point_in_ray(const double3& ray_start, const double3& ray_end,
                         const double t);
                         
#endif // VOXEL_TRAVERSAL_H
