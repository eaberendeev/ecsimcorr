// Author: Evgeny Berendeev
// Budker Institute of Nuclear Physics of Siberian Branch Russian Academy of Sciences
// beren@inp.nsk.su
// (c) 2022, for licensing details see the LICENSE file

#pragma once

#ifndef UTIL_H
#define UTIL_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <unordered_map>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <algorithm>

enum Axis { X = 0, Y, Z };

enum CommTags { PARTICLES = 0, FIELDS };

// Define basic types
#define MAJOR Eigen::RowMajor 
typedef Eigen::SparseMatrix<double, MAJOR> Operator;
typedef Eigen::VectorXd Field;
typedef Eigen::Vector2i Vector2i;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Triplet<double> Trip;
typedef std::unordered_map<int,double> IndexMap;

//general indexing routine (row major)
inline constexpr int ind(int x, int y, int z, int c, int Nx, int Ny, int Nz, int Nc)
{
	return(c + Nc*(z + Nz*(y + Ny*x)));
}

enum status{

};
#endif
