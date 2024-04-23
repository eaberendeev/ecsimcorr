#ifndef SOLVERSLE_H_
#define SOLVERSLE_H_

#include "util.h"
#include "containers.h"
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/IterativeLinearSolvers>

#define SOLVER_TOLERANCE 1.e-6
// System of linear equations solvers
typedef Eigen::GMRES<Eigen::SparseMatrix<double, MAJOR> > gmres;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double, MAJOR> > bicgstab;

void solve_SLE(const Operator &A, const Field& rhs, Field& x, const Field& x0);
#endif 	
