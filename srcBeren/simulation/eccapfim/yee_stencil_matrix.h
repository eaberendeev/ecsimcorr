#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>
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

    // Preconditioners of the GMRES solve:
    // - BlockJacobi: inverse of the 3x3 blocks of the three E components of a node;
    // - Circulant: exact inverse, by FFT, of the operator averaged over all
    //   translations (T. Chan's circulant approximation): the Maxwell part is
    //   translation invariant, the plasma part is replaced by its mean. Exact
    //   for a uniform plasma; needs a fully periodic grid without constraints;
    // - CirculantBlockJacobi: symmetric two-level product BJ, Circulant, BJ,
    //   which adds the local density variation.
    enum class Preconditioner { BlockJacobi, Circulant, CirculantBlockJacobi };

    // `periodic`: every axis is periodic and every unknown is free.
    void init(const Eigen::Vector3i& first_node, const Eigen::Vector3i& nodes, const DofMap& dof,
              const Vector& free_rows, const Eigen::SparseMatrix<double>& constant_part, bool periodic);

    bool supports_circulant() const {
        return periodic_;
    }

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

    // The same without atomic update, for callers that keep rows disjoint.
    void add_unsynchronized(int row, int c1, int c2, int ox, int oy, int oz, double value) {
        const auto& box = boxes_[c1][c2];
        const size_t slot = static_cast<size_t>(row) * kSlots + box.base +
                            ((ox - box.lo[0]) * box.width[1] + (oy - box.lo[1])) * box.width[2] + (oz - box.lo[2]);
        if (enabled_[slot])
            values_[slot] += value;
    }

    int rows() const {
        return rows_;
    }

    void multiply(const Vector& x, Vector& y) const;
    // Prepares the preconditioner from the current values.
    void factor(Preconditioner preconditioner);
    void apply_preconditioner(const Vector& x, Vector& y) const;
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
    void factor_block_jacobi();
    void apply_block_jacobi(const Vector& x, Vector& y) const;
    void factor_circulant();
    void apply_circulant(const Vector& x, Vector& y) const;
    Preconditioner preconditioner_ = Preconditioner::BlockJacobi;
    bool periodic_ = false;
    int nodes_[3] = {0, 0, 0};
    std::vector<Eigen::Matrix3cd> symbol_inverse_;   // per Fourier mode
};
