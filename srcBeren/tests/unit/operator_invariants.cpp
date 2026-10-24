#include <omp.h>

#include <random>

#include "Damping.h"
#include "sparse.h"
#include "test_support.h"
#include "thread_partitioned_matrix.h"

void sparse_algebra() {
    for (int scenario = 0; scenario < 8; ++scenario) {
        const int rows = scenario == 0 ? 0 : 19, cols = 23;
        Operator a(rows, cols), b(rows, cols);
        std::vector<Eigen::Triplet<double>> ta, tb;
        if (scenario > 1)
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c) {
                    if ((r * 7 + c * 3 + scenario) % 11 == 0)
                        ta.emplace_back(r, c, 0.25 * (r - c));
                    if ((r + c + scenario) % 7 == 0)
                        tb.emplace_back(r, c, -0.25 * (r - c));
                }
        a.setFromTriplets(ta.begin(), ta.end());
        b.setFromTriplets(tb.begin(), tb.end());
        const Operator expected = a + b, actual = parallelSparseSum(a, b);
        if (rows > 0)
            check::near((expected - actual).norm(), 0, 1e-13, "sparse union: empty rows, overlap, cancellation");
        else
            check::require(actual.rows() == 0 && actual.cols() == cols && actual.nonZeros() == 0,
                           "empty sparse matrix shape/portrait changed");
        Field3d x(cols), result(rows);
        for (int i = 0; i < cols; ++i) x[i] = std::cos(i * 0.7);
        spmv(actual, x, result);
        check::near((result.data() - expected * x.data()).norm(), 0, 1e-12, "CSR multiply vs independent Eigen path");
        if (rows > 0) {
            ThreadPartitionedSparseMatrix<double> partitioned(actual);
            spmv(partitioned, x, result);
            check::near((result.data() - expected * x.data()).norm(), 0, 1e-12, "partitioned CSR multiply");
        }
    }
}

void reductions() {
    for (int size : {0, 1, 31, 16383, 32769}) {
        Eigen::VectorXd x(size), y(size);
        long double dot = 0, square = 0;
        for (int i = 0; i < size; ++i) {
            x[i] = std::sin(i * 0.3);
            y[i] = std::cos(i * 0.7);
            dot += static_cast<long double>(x[i]) * y[i];
            square += static_cast<long double>(x[i]) * x[i];
        }
        for (int threads : {1, 2, 4}) {
            omp_set_num_threads(threads);
            check::near(blas::dot(x, y), double(dot), 1e-9, "SIMD/OpenMP dot");
            check::near(blas::squaredNorm(x), double(square), 1e-8, "SIMD/OpenMP norm");
            check::near(blas::normalizedDot(x, y), square ? double(dot / square) : 0, 1e-12,
                        "SIMD/OpenMP normalized dot");
            const auto first = blas::dot(x, y);
            check::near(blas::dot(x, y), first, 0, "fixed-thread reductions must repeat exactly");
        }
    }
    omp_set_num_threads(2);
}

void maxwell_operators() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(0.7, 1.1, 0.9));
    Mesh mesh;
    BoundaryConditionHandler bc;
    const int n = 3 * domain.total_size();
    Operator curl_e(n, n), curl_b(n, n), div(domain.total_size(), n);
    mesh.stencil_curlE(curl_e, domain, bc);
    mesh.stencil_curlB(curl_b, domain, bc);
    mesh.stencil_divE(div, domain, bc);
    Field3d e(domain.size(), 3), b(domain.size(), 3);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                for (int d = 0; d < 3; ++d) {
                    e(i, j, k, d) = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d)[(d + 1) % 3];
                    b(i, j, k, d) = domain.get_node_position(i, j, k, FieldType::MAGNETIC, d)[(d + 1) % 3];
                }
    const Eigen::VectorXd ce = curl_e * e.data(), cb = curl_b * b.data();
    for (int i = 2; i < domain.size().x() - 2; ++i)
        for (int j = 2; j < domain.size().y() - 2; ++j)
            for (int k = 2; k < domain.size().z() - 2; ++k)
                for (int d = 0; d < 3; ++d) {
                    check::near(ce[domain.vind(i, j, k, d)], -1, 1e-13, "curl E on an affine field");
                    check::near(cb[domain.vind(i, j, k, d)], -1, 1e-13, "curl B on an affine field");
                }
    std::mt19937 engine(17);
    std::uniform_real_distribution<double> uniform(-1, 1);
    for (int i = 0; i < e.capacity(); ++i) e[i] = uniform(engine);
    const Eigen::VectorXd identity = div * (curl_b * e.data());
    for (int i = 2; i < domain.size().x() - 2; ++i)
        for (int j = 2; j < domain.size().y() - 2; ++j)
            for (int k = 2; k < domain.size().z() - 2; ++k)
                check::near(identity[domain.sind(i, j, k)], 0, 1e-12, "discrete div(curl B) must vanish");
}

void damping_balance() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(1, 1, 1));
    for (const std::string type : {"Rectangle", "CircleXY"}) {
        Field3d e(domain.size(), 3), b(domain.size(), 3);
        e.data().setConstant(0.7);
        b.data().setConstant(-0.3);
        const double before = 0.5 * (e.squared() + b.squared());
        nlohmann::json cfg = {{"DampingType", type},
                              {"DampingCoeff", 0.8},
                              {"DampCellsX_glob", {2, 1}},
                              {"DampCellsY_glob", {1, 2}},
                              {"DampCellsZ_glob", {2, 2}}};
        const double removed = damping_fields(e, b, domain, cfg);
        check::near(removed, before - 0.5 * (e.squared() + b.squared()), 1e-10,
                    "damping must account for all removed energy");
        check::require(removed > 0, "nonempty damping layer did nothing");
    }
}
int main() {
    return check::run([] {
        sparse_algebra();
        reductions();
        maxwell_operators();
        damping_balance();
    });
}
