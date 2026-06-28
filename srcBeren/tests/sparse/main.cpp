#include <omp.h>

#include <Eigen/Sparse>
#include <chrono>
#include <iostream>
#include <random>

#include "sparse.h"
#include "timer.h"

using SparseMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;

SparseMat genSparseMatrix(int64_t rows, int64_t cols, int64_t nnz, int64_t seed) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> rowDistr(0, rows);
    std::uniform_int_distribution<int> colDistr(0, cols);
    std::uniform_real_distribution<double> valueDistr(-1.0, 1.0);

    std::vector<Eigen::Triplet<int>> triplets(nnz);

    for (int64_t i = 0; i < nnz; ++i) {
        const Eigen::Triplet<int> triplet(rowDistr(gen), colDistr(gen), valueDistr(gen));
        triplets[i] = triplet;
    }

    SparseMat res(rows, cols);
    res.setFromTriplets(triplets.begin(), triplets.end());

    return res;
}

int main() {
    std::cout << "test of optimized parallel and reference summation (from Eigen) of sparse matrices" << std::endl;
    std::cout << "Number of omp threads: " << omp_get_max_threads() << std::endl;

    constexpr int64_t rows = 1'000'000;
    constexpr int64_t cols = 2'000'000;
    constexpr int64_t nnzA = rows * 81;
    constexpr int64_t nnzB = rows * 150;
    constexpr int64_t nnzC = rows * 15;

    std::cout << "test of A + BC, with" << std::endl;
    std::cout << "rows: " << rows << std::endl;
    std::cout << "cols: " << cols << std::endl;
    std::cout << "nnzA: " << nnzA << std::endl;
    std::cout << "nnzBC in range [" << std::min(nnzB, nnzC) << "; " << nnzB + nnzC << "]" << std::endl;

    const SparseMat A = genSparseMatrix(rows, cols, nnzA, 32);
    const SparseMat B = genSparseMatrix(rows, cols, nnzB, 33);
    const SparseMat C = genSparseMatrix(rows, cols, nnzC, 34);

    // We test 2 matrices, which have intersection elements, so test performance of A - BC
    const SparseMat BC = B + C;

    const double time1 = omp_get_wtime();
    const SparseMat refSum = A + BC;
    const double time2 = omp_get_wtime();
    const SparseMat testSum = parallelSparseSum(A, BC);
    const double time3 = omp_get_wtime();

    const SparseMat diff = refSum - testSum;

    const double err = diff.norm();

    std::cout << "ref time:  " << time2 - time1 << "[s]" << std::endl;
    std::cout << "test time: " << time3 - time2 << "[s]" << std::endl;
    std::cout << "||ref - test||: " << err << " (expected 0)" << std::endl;

    return (err == 0.0) ? 0 : 1;
}
