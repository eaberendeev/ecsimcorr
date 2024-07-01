#include <cfloat>
#include <iostream>
#include <vector>
#include "Vec.h"

#include "voxel_traversal.h"

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
                                             const double3& ray_end, const double bin_size) {
  std::vector<int3> visited_voxels;

  // This id of the first/current voxel hit by the ray.
  // Using floor (round down) is actually very important,
  // the implicit int-casting will round up for negative numbers.
  int3 current_voxel(std::floor(ray_start(0) / bin_size),
                     std::floor(ray_start(1) / bin_size),
                     std::floor(ray_start(2) / bin_size));

  // The id of the last voxel hit by the ray.
  // TODO: what happens if the end point is on a border?
  int3 last_voxel(std::floor(ray_end(0) / bin_size),
                  std::floor(ray_end(1) / bin_size),
                  std::floor(ray_end(2) / bin_size));

  if (current_voxel == last_voxel) {
    visited_voxels.push_back(current_voxel);
    return visited_voxels;
  }

      // Compute normalized ray direction.
      double3 ray = ray_end - ray_start;
  // ray.normalize();

  // In which direction the voxel ids are incremented.
  double stepX = (ray(0) >= 0) ? 1 : -1; // correct
  double stepY = (ray(1) >= 0) ? 1 : -1; // correct
  double stepZ = (ray(2) >= 0) ? 1 : -1; // correct

  // Distance along the ray to the next voxel border from the current position
  // (tMaxX, tMaxY, tMaxZ).
  double next_voxel_boundary_x =
      (current_voxel(0) + stepX) * bin_size; // correct
  double next_voxel_boundary_y =
      (current_voxel(1) + stepY) * bin_size; // correct
  double next_voxel_boundary_z =
      (current_voxel(2) + stepZ) * bin_size; // correct

  next_voxel_boundary_z += (stepZ < 0) ? bin_size : 0;
  next_voxel_boundary_y += (stepY < 0) ? bin_size : 0;
  next_voxel_boundary_x += (stepX < 0) ? bin_size : 0;
  // tMaxX, tMaxY, tMaxZ -- distance until next intersection with voxel-border
  // the value of t at which the ray crosses the first vertical voxel boundary
  double tMaxX = (ray(0) != 0) ? (next_voxel_boundary_x - ray_start(0)) / ray(0)
                               : DBL_MAX; //
  double tMaxY = (ray(1) != 0) ? (next_voxel_boundary_y - ray_start(1)) / ray(1)
                               : DBL_MAX; //
  double tMaxZ = (ray(2) != 0) ? (next_voxel_boundary_z - ray_start(2)) / ray(2)
                               : DBL_MAX; //

  // tDeltaX, tDeltaY, tDeltaZ --
  // how far along the ray we must move for the horizontal component to equal
  // the width of a voxel the direction in which we traverse the grid can only
  // be FLT_MAX if we never go in that direction
  double tDeltaX = (ray(0) != 0) ? bin_size / ray(0) * stepX : DBL_MAX;
  double tDeltaY = (ray(1) != 0) ? bin_size / ray(1) * stepY : DBL_MAX;
  double tDeltaZ = (ray(2) != 0) ? bin_size / ray(2) * stepZ : DBL_MAX;

  int3 diff(0, 0, 0);
  bool neg_ray = false;
  for(int dim = 0; dim < 3; dim++){
    if(current_voxel(dim) != last_voxel(dim) && ray(dim) < 0){
      diff(dim)--;
      neg_ray = true;
    }
  }

  visited_voxels.push_back(current_voxel);


  while (last_voxel != current_voxel) {
    if (tMaxX < tMaxY) {
      if (tMaxX < tMaxZ) {
        current_voxel(0) += stepX;
        tMaxX += tDeltaX;
      } else {
        current_voxel(2) += stepZ;
        tMaxZ += tDeltaZ;
      }
    } else {
      if (tMaxY < tMaxZ) {
        current_voxel(1) += stepY;
        tMaxY += tDeltaY;
      } else {
        current_voxel(2) += stepZ;
        tMaxZ += tDeltaZ;
      }
    }
    visited_voxels.push_back(current_voxel);
  }
  return visited_voxels;
}

double find_ray_voxel_intersection_parameter(const double3& ray_start,
                                                  const double3& ray_end,
                                                  const int3& current_voxel,
                                                  const int3& next_voxel,
                                                  double bin_size) {
    const int3 diff = next_voxel - current_voxel;
    int ax, dir;
    if (diff(0) != 0) {
        ax = 0;
    } else if (diff(1) != 0) {
        ax = 1;
    } else if (diff(2) != 0) {
        ax = 2;
    } else {
        return 1;
    }
    dir = diff(ax) > 0 ? 1 : 0;

    double t = ((current_voxel(ax) + dir) * bin_size - ray_start(ax)) /
               (ray_end(ax) - ray_start(ax));
    return t;
}

double3 get_point_in_ray(const double3& ray_start, const double3& ray_end,
                          const double t) {
  double3 point;
  for (int i = 0; i < 3; i++) {
    point(i) = ray_start(i) + t * (ray_end(i) - ray_start(i));
  }
  return point;
}


// int main(int, char **) {
//   double3 ray_start(1.4, 0.6, 0.1);
//   double3 ray_end(0.2, 0.5, 0.1);
//   double bin_size = 0.1;
//   std::cout << "Voxel size: " << bin_size << std::endl;
//   std::cout << "Starting position: " << ray_start << std::endl;
//   std::cout << "Ending position: " << ray_end << std::endl;
//   std::cout << "Voxel ID's from start to end:" << std::endl;
//   std::vector<int3> ids = voxel_traversal(ray_start, ray_end, bin_size);

//   for (auto &i : ids) {
//     std::cout << "> " << i << std::endl;
//   }
//   std::cout << "Total number of traversed voxels: " << ids.size() << std::endl;
//   int3 current_voxel = ids(0);
//   for (int i =1; i < ids.size(); i++) {
//     double t = find_t(ray_start, ray_end, current_voxel, ids(i), bin_size);
//     double3 point = get_point(ray_start, ray_end, t);
//     std::cout << "t: " << t << " "<< point <<std::endl;
//     current_voxel = ids(i);
//   }
//     return 0;
// }
