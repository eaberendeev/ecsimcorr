// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef UTIL_H
#define UTIL_H
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace Axis{
	enum { X = 0, Y, Z };
};
namespace Dim{
	enum {
    X = 0,
    Y,
    Z,
    COUNT   // Avoid hardcoded value
};
};
constexpr auto MAX_DIM = Dim::COUNT;

const int SHAPE_SIZE = 2;
// ghost cells for each side
const int GHOST_CELLS = 1;
/**
 * nodes = cells + 1
 * GHOST_NODES = 2 * GHOST_CELLS + 1
*/
const int GHOST_NODES = 3;

enum CommTags { PARTICLES = 0, FIELDS };
enum BoundType { PERIODIC = 0, OPEN, NEIGHBOUR };

// Define basic types
#define MAJOR Eigen::RowMajor
typedef Eigen::SparseMatrix<double, MAJOR> Operator;
typedef Eigen::VectorXd Field;
typedef Eigen::Vector2i Vector2i;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Triplet<double> Trip;
typedef std::unordered_map<int, double> IndexMap;

// general indexing routine (row major)
inline constexpr int ind(int x, int y, int z, int c, int Nx, int Ny, int Nz,
                         int Nc) {
    return (c + Nc * (z + Nz * (y + Ny * x + Nx * 0)));
}

enum status {

};
#endif
