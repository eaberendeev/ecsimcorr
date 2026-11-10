#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <vector>

// Linear operator on the reduced E unknowns with the fixed coupling pattern
// of one particle: E_c uses S1 along c and S2 across c, so the product of the
// weights of E_c1 and E_c2 couples nodes whose offset per axis a lies in
//   [-1, 1]  a = c1 = c2,   [-1, 2]  a = c1 != c2,
//   [-2, 1]  a = c2 != c1,  [-2, 2]  otherwise.
// That is 75 + 80 + 80 = 235 slots per row. The Yee curl-curl operator fits
// into the same pattern, so A = I + dt^2/4 curl_B curl_E + dt/2 M_p is one
// array of values with a precomputed column per slot. Slots that leave the
// unknowns (Dirichlet columns, nonperiodic faces) point to their own row and
// stay zero.
class YeeStencilMatrix {
   public:
    static constexpr int kSlots = 235;
    using Vector = Eigen::VectorXd;
    // Reduced index of component c at storage node (i, j, k), or -1.
    using DofMap = std::function<int(int i, int j, int k, int c)>;

    void init(const Eigen::Vector3i& first_node, const Eigen::Vector3i& nodes, const DofMap& dof,
              const Vector& free_rows, const Eigen::SparseMatrix<double>& constant_part);

    // values = constant part (the Maxwell operator), ready for particle terms.
    void reset();

    // Adds `value` to `row` (component c1) and the column of component c2 at
    // the node offset (ox, oy, oz) from the row node, the unwrapped offset of
    // the two weights. Disabled slots (constrained columns) are skipped.
    void add(int row, int c1, int c2, int ox, int oy, int oz, double value) {
        const auto& box = boxes_[c1][c2];
        const size_t slot = static_cast<size_t>(row) * kSlots + box.base +
                            ((ox - box.lo[0]) * box.width[1] + (oy - box.lo[1])) * box.width[2] + (oz - box.lo[2]);
        if (!enabled_[slot])
            return;
#pragma omp atomic update
        values_[slot] += value;
    }

    int rows() const {
        return rows_;
    }

    void multiply(const Vector& x, Vector& y) const;
    // Applies the inverse of the 3x3 node blocks of the current values.
    void factor_block_jacobi();
    void apply_block_jacobi(const Vector& x, Vector& y) const;
    // Right-preconditioned restarted GMRES; returns false if not converged.
    bool solve_gmres(const Vector& rhs, Vector& x, double tolerance, int restart, int max_iterations,
                     int& iterations) const;
    Eigen::SparseMatrix<double> to_sparse() const;

   private:
    struct Box {
        int base = 0;
        int lo[3] = {0, 0, 0};
        int width[3] = {0, 0, 0};
    };

    Box boxes_[3][3];
    int rows_ = 0;
    std::vector<int> columns_;
    std::vector<double> values_, constant_values_;
    std::vector<std::vector<std::pair<int, int>>> block_slots_;   // per row: (slot, c2) inside the node block
    std::vector<Eigen::Matrix3d> block_inverse_;
    std::vector<char> free_, enabled_;
};
