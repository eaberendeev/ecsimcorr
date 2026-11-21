#pragma once

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCore>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace eccapfim {

// Block-Jacobi preconditioner for a field matrix laid out as three consecutive
// components per node. It keeps only the 3x3 diagonal node blocks, which is the
// cheapest separable approximation to the Maxwell operator.
template <typename MatrixType_>
class BlockJacobi {
   public:
    using MatrixType = MatrixType_;
    using Scalar = typename MatrixType::Scalar;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using StorageIndex = typename Vector::StorageIndex;
    enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic };

    BlockJacobi() : m_isInitialized(false) {}

    template <typename MatType>
    explicit BlockJacobi(const MatType& mat) : m_isInitialized(false) {
        compute(mat);
    }

    Eigen::Index rows() const {
        return 3 * static_cast<Eigen::Index>(m_inverse.size());
    }
    Eigen::Index cols() const {
        return rows();
    }

    template <typename MatType>
    BlockJacobi& compute(const MatType& mat) {
        if (mat.rows() != mat.cols() || mat.rows() % 3 != 0)
            throw std::invalid_argument("block-Jacobi requires a square matrix with three components per node");
        const int blocks = static_cast<int>(mat.rows() / 3);
        std::vector<Eigen::Matrix<Scalar, 3, 3>> local(blocks, Eigen::Matrix<Scalar, 3, 3>::Zero());
        for (int outer = 0; outer < mat.outerSize(); ++outer)
            for (typename MatType::InnerIterator it(mat, outer); it; ++it) {
                const int row = static_cast<int>(it.row()), col = static_cast<int>(it.col());
                if (row / 3 == col / 3)
                    local[row / 3](row % 3, col % 3) = it.value();
            }
        m_inverse.resize(blocks);
        for (int b = 0; b < blocks; ++b) {
            Eigen::FullPivLU<Eigen::Matrix<Scalar, 3, 3>> lu(local[b]);
            if (lu.isInvertible()) {
                m_inverse[b] = lu.inverse();
            } else {
                m_inverse[b].setZero();
                for (int i = 0; i < 3; ++i)
                    if (local[b](i, i) != Scalar(0))
                        m_inverse[b](i, i) = Scalar(1) / local[b](i, i);
            }
        }
        m_isInitialized = true;
        return *this;
    }

    template <typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const {
        x.resize(b.rows());
        if (b.rows() != rows())
            throw std::invalid_argument("block-Jacobi rhs size mismatch");
        for (Eigen::Index i = 0; i < b.rows(); i += 3) {
            Eigen::Matrix<Scalar, 3, 1> rhs(b[i], b[i + 1], b[i + 2]);
            const Eigen::Matrix<Scalar, 3, 1> solution = m_inverse[static_cast<size_t>(i / 3)] * rhs;
            x[i] = solution[0];
            x[i + 1] = solution[1];
            x[i + 2] = solution[2];
        }
    }

    template <typename Rhs>
    inline const Eigen::Solve<BlockJacobi, Rhs> solve(const Eigen::MatrixBase<Rhs>& b) const {
        if (!m_isInitialized)
            throw std::logic_error("block-Jacobi preconditioner is not initialized");
        return Eigen::Solve<BlockJacobi, Rhs>(*this, b.derived());
    }

    Eigen::ComputationInfo info() const {
        return Eigen::Success;
    }

   private:
    std::vector<Eigen::Matrix<Scalar, 3, 3>> m_inverse;
    bool m_isInitialized;
};

}   // namespace eccapfim
