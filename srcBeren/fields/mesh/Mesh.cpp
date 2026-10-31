#include "Mesh.h"

#include "coupling/shape.h"
#include "grid/domain.h"
#include "coupling/gather.h"
#include "row_block.h"
#include "solverSLE.h"
#include "thread_partitioned_matrix.h"
#include "timer.h"

void Mesh::init(const Domain& domain, double dt, BoundaryConditionHandler& bc_handler) {
    RECORD_TIMER;

    Lmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Lmat2.resize(domain.total_size() * 3, domain.total_size() * 3);
    Mmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    Imat.resize(domain.total_size() * 3, domain.total_size() * 3);
    IMmat.resize(domain.total_size() * 3, domain.total_size() * 3);
    divE.resize(domain.total_size(), domain.total_size() * 3);

    xCellSize = domain.cell_size().x();
    yCellSize = domain.cell_size().y();
    zCellSize = domain.cell_size().z();
    xSize = domain.size().x();
    ySize = domain.size().y();
    zSize = domain.size().z();

    Operator curlBtmp;
    Operator curlEtmp;
    curlBtmp.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlEtmp.resize(domain.total_size() * 3, domain.total_size() * 3);

    stencil_Imat(Imat, domain);
    stencil_curlE(curlEtmp, domain, bc_handler);
    stencil_curlB(curlBtmp, domain, bc_handler);

    curlE = ThreadPartitionedSparseMatrix<double>(curlEtmp);
    curlB = ThreadPartitionedSparseMatrix<double>(curlBtmp);

    stencil_divE(divE, domain, bc_handler);

    Mmat = -0.25 * dt * dt * curlBtmp * curlEtmp;
    IMmat = Imat - Mmat;
    IMmat.makeCompressed();
    // Дирихле-строки для предписанных полей (er0 mode "solv") закрепляются
    // один раз, т.к. IMmat статическая; er0-строки Lmat2 занулены row-only
    // фильтром, поэтому A = IMmat + Lmat2 сохраняет пиннинг, а solve_electric_field_from_deposited_current
    // использует IMmat напрямую.
    impose_dirichlet_rows(IMmat, bc_handler.electric_dirichlet(domain));
}

void Mesh::print_operator(const Operator& oper) {
    for (int k = 0; k < oper.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(oper, k); it; ++it) {
            std::cout << pos_vind(it.row(), 0) << " " << pos_vind(it.row(), 1) << " " << pos_vind(it.row(), 2) << " "
                      << pos_vind(it.row(), 3) << " " << pos_vind(it.col(), 0) << " " << pos_vind(it.col(), 1) << " "
                      << pos_vind(it.col(), 2) << " " << pos_vind(it.col(), 3) << " " << it.value() << "\n";
        }
    }
}

void Mesh::prepare() {
}

void Mesh::fdtd_explicit(Field3d& E, Field3d& B, const Field3d& J, const double dt) {
    RECORD_TIMER;
    E += 0.5 * dt * (curlB * B) - 0.5 * dt * J;
    B -= 0.5 * dt * (curlE * E);
}

void Mesh::computeB(const Field3d& fieldE, const Field3d& fieldEn, Field3d& fieldB, double dt) {
    RECORD_TIMER;
    fieldB -= (0.5 * dt) * (curlE * (fieldE + fieldEn));
}

void Mesh::compute_fieldB(Field3d& Bn, const Field3d& B, const Field3d& E, const Field3d& En, double dt) {
    RECORD_TIMER;
    Bn = B - (0.5 * dt) * (curlE * (E + En));
}
