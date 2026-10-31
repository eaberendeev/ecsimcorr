#include "maxwell_operators.h"

#include <utility>

#include "timer.h"

namespace {

void assemble_curlB(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    RECORD_TIMER;

    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    auto addRowIfInside = [&](int i, int j, int k, int comp) -> bool {
        return domain.is_inside_node_periodic(i, j, k, FieldType::ELECTRIC, comp, is_periodic);
    };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int im = bc_handler.wrap_index(i - 1, 0, domain);
                const int jm = bc_handler.wrap_index(j - 1, 1, domain);
                const int km = bc_handler.wrap_index(k - 1, 2, domain);

                const int vindx = domain.vind(i, j, k, 0);   // Ex
                const int vindy = domain.vind(i, j, k, 1);   // Ey
                const int vindz = domain.vind(i, j, k, 2);   // Ez

                // (x)[i+1/2,j,k]
                if (addRowIfInside(i, j, k, 0)) {
                    // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                    double val = 1.0 / dy;
                    trips.emplace_back(vindx, domain.vind(i, j, k, 2), val);
                    trips.emplace_back(vindx, domain.vind(i, jm, k, 2), -val);
                    // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                    val = -1.0 / dz;
                    trips.emplace_back(vindx, domain.vind(i, j, k, 1), val);
                    trips.emplace_back(vindx, domain.vind(i, j, km, 1), -val);
                }
                // (y)[i,j+1/2,k]
                if (addRowIfInside(i, j, k, 1)) {
                    // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                    double val = 1.0 / dz;
                    trips.emplace_back(vindy, domain.vind(i, j, k, 0), val);
                    trips.emplace_back(vindy, domain.vind(i, j, km, 0), -val);
                    // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                    val = -1.0 / dx;
                    trips.emplace_back(vindy, domain.vind(i, j, k, 2), val);
                    trips.emplace_back(vindy, domain.vind(im, j, k, 2), -val);
                }
                // (z)[i,j,k+1/2]
                if (addRowIfInside(i, j, k, 2)) {
                    // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                    double val = 1.0 / dx;
                    trips.emplace_back(vindz, domain.vind(i, j, k, 1), val);
                    trips.emplace_back(vindz, domain.vind(im, j, k, 1), -val);
                    // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                    val = -1.0 / dy;
                    trips.emplace_back(vindz, domain.vind(i, j, k, 0), val);
                    trips.emplace_back(vindz, domain.vind(i, jm, k, 0), -val);
                }
                //  bc_handler.modify_curlB_stencil(i, j, k, trips, domain);
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void assemble_curlE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    RECORD_TIMER;

    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    // Проверяет, принадлежит ли магнитный узел области (с учётом периодичности)
    auto addRowIfInside = [&](int i, int j, int k, int comp) -> bool {
        return domain.is_inside_node_periodic(i, j, k, FieldType::MAGNETIC, comp, is_periodic);
    };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int ip = bc_handler.wrap_index(i + 1, 0, domain);
                const int jp = bc_handler.wrap_index(j + 1, 1, domain);
                const int kp = bc_handler.wrap_index(k + 1, 2, domain);

                const int vindx = domain.vind(i, j, k, 0);   // Bx
                const int vindy = domain.vind(i, j, k, 1);   // By
                const int vindz = domain.vind(i, j, k, 2);   // Bz

                // (x)[i,j+1/2,k+1/2]
                if (addRowIfInside(i, j, k, 0)) {
                    // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
                    double val = 1.0 / dy;
                    trips.emplace_back(vindx, domain.vind(i, jp, k, 2), val);
                    trips.emplace_back(vindx, domain.vind(i, j, k, 2), -val);
                    // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                    val = -1.0 / dz;
                    trips.emplace_back(vindx, domain.vind(i, j, kp, 1), val);
                    trips.emplace_back(vindx, domain.vind(i, j, k, 1), -val);
                }

                // (y)[i+1/2,j,k+1/2]
                if (addRowIfInside(i, j, k, 1)) {
                    // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                    double val = 1.0 / dz;
                    trips.emplace_back(vindy, domain.vind(i, j, kp, 0), val);
                    trips.emplace_back(vindy, domain.vind(i, j, k, 0), -val);
                    // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                    val = -1.0 / dx;
                    trips.emplace_back(vindy, domain.vind(ip, j, k, 2), val);
                    trips.emplace_back(vindy, domain.vind(i, j, k, 2), -val);
                }

                // (z)[i+1/2,j+1/2,k]
                if (addRowIfInside(i, j, k, 2)) {
                    // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                    double val = 1.0 / dx;
                    trips.emplace_back(vindz, domain.vind(ip, j, k, 1), val);
                    trips.emplace_back(vindz, domain.vind(i, j, k, 1), -val);
                    // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                    val = -1.0 / dy;
                    trips.emplace_back(vindz, domain.vind(i, jp, k, 0), val);
                    trips.emplace_back(vindz, domain.vind(i, j, k, 0), -val);
                }

                //  bc_handler.modify_curlE_stencil(i, j, k, trips, domain);
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void assemble_divE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    RECORD_TIMER;

    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 6;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                if (!domain.is_inside_node_periodic(i, j, k, FieldType::DENSITY, 0, is_periodic))
                    continue;

                const int im = bc_handler.wrap_index(i - 1, 0, domain);
                const int jm = bc_handler.wrap_index(j - 1, 1, domain);
                const int km = bc_handler.wrap_index(k - 1, 2, domain);

                const int sindx = domain.sind(i, j, k);

                // [i,j,k]
                // ( Ex[i+1/2,j,k] - Ex[i-1/2,j,k] ) / dx
                double val = 1.0 / dx;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 0), val));
                trips.push_back(Trip(sindx, domain.vind(im, j, k, 0), -val));
                // ( Ey[i,j+1/2,k] - Ey[i,j-1/2,k] ) / dy
                val = 1.0 / dy;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 1), val));
                trips.push_back(Trip(sindx, domain.vind(i, jm, k, 1), -val));
                // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
                val = 1.0 / dz;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 2), val));
                trips.push_back(Trip(sindx, domain.vind(i, j, km, 2), -val));
                // TODO: add bc_handler modify stencil
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

}   // namespace

MaxwellStencils assemble_maxwell_stencils(const Domain& domain, BoundaryConditionHandler& boundaries) {
    RECORD_TIMER;
    const int n = 3 * domain.total_size();
    MaxwellStencils stencils{Operator(n, n), Operator(n, n), Operator(domain.total_size(), n)};
    assemble_curlE(stencils.curlE, domain, boundaries);
    assemble_curlB(stencils.curlB, domain, boundaries);
    assemble_divE(stencils.divE, domain, boundaries);
    return stencils;
}

void MaxwellOperators::init(MaxwellStencils&& stencils) {
    RECORD_TIMER;
    curlE = ThreadPartitionedSparseMatrix<double>(stencils.curlE);
    curlB = ThreadPartitionedSparseMatrix<double>(stencils.curlB);
    divE = ThreadPartitionedSparseMatrix<double>(stencils.divE);
}

void MaxwellOperators::advance_magnetic_field(Field3d& Bn, const Field3d& B, const Field3d& E, const Field3d& En,
                                              double dt) const {
    RECORD_TIMER;
    Bn = B - (0.5 * dt) * (curlE * (E + En));
}
