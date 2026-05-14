#pragma once

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <unordered_map>
#include <unsupported/Eigen/IterativeSolvers>

enum Axis : int {
    X = 0,
    Y = 1,
    Z = 2,
    C = 3,
};

enum class ShapeType { NGP, Linear, Quadratic };

enum Preconditioner { NONE, DIAGONAL };

#define MAJOR Eigen::RowMajor
typedef Eigen::SparseMatrix<double, MAJOR> Operator;
typedef Eigen::GMRES<Eigen::SparseMatrix<double, MAJOR>> gmres;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double, MAJOR>> bicgstab;
typedef Eigen::VectorXd Field;
typedef Eigen::Vector2i Vector2i;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Triplet<double> Trip;
typedef std::unordered_map<int, double> IndexMap;
