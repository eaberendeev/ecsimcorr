#include <memory>
#include <source_location>

#include "Mesh.h"
#include "coupling/shape.h"
#include "config.h"
#include "env_options.h"
#include "grid/domain.h"
#include "log_macros.h"
#include "memory.h"
#include "timer.h"
#include "util.h"

// matrix ColMajor
// (0,0) (0,1) (0,2)
// (1,0) (1,1) (1,2)
// (2,0) (2,1) (2,2)

// Ex(i+/2,j,k), Ey(i,j+1/2,k), Ez(i,j,k+1/2)
// Bx(i,j+1/2,k+1/2), By(i+1/2,j,k+1/2), Bz(i+/2,j+1/2,k)

Mesh::Mesh() = default;
Mesh::~Mesh() = default;

void Mesh::stencil_Imat(Operator& mat, const Domain& domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z();
    trips.reserve(totalSize);

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 0), vind(i, j, k, 0), 1.0));

                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 1), vind(i, j, k, 1), 1.0));

                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 2), vind(i, j, k, 2), 1.0));
            }
        }
    }
    LOG_STEP("Stencil_Imat done. size: " << size << ". Elements: " << 3 * totalSize << " " << Imat.rows() << " "
                                         << Imat.cols() << "\n");
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlB(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
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

void Mesh::stencil_curlE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
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

void Mesh::stencil_divE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
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

void Mesh::stencil_smooth_1d(Operator& mat, const Domain& domain, int dim) {
    RECORD_TIMER;

    constexpr int COMPONENTS = 3;
    constexpr int STENCIL = 3;
    const auto size = domain.size();
    const int nx = size.x();
    const int ny = size.y();
    const int nz = size.z();

    const int totalCells = nx * ny * nz;
    // checker
    mat.resize(domain.total_size() * 3, domain.total_size() * 3);

    if (dim != static_cast<int>(Axis::X) && dim != static_cast<int>(Axis::Y) && dim != static_cast<int>(Axis::Z)) {
        std::cout << "Error: Invalid dimension." << std::endl;
        return;
    }
    std::vector<Trip> trips;
    trips.reserve(static_cast<size_t>(totalCells * COMPONENTS * STENCIL));

    auto addTriplet = [](std::vector<Trip>& trips, const Domain& domain, int row, int col, double val) {
        bool onArea = domain.is_inside_node(row, FieldType::ELECTRIC);
        if (onArea) {
            trips.push_back(Trip(row, col, val));
        }
    };

    // auto imod = [](int a, int m) {
    //     // % m - 1
    //     if (a > m - 1) return 3;
    //     if (a < 0) return m - 4;
    //     return a;
    // };

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                int im = i - 1, ip = i + 1;
                int jm = j - 1, jp = j + 1;
                int km = k - 1, kp = k + 1;

                // if (dim == static_cast<int>(Axis::X) && domain.is_periodic_bound(X)) {
                //     im = imod(i - 1, nx);
                //     ip = imod(i + 1, nx);
                // } else if (dim == static_cast<int>(Axis::Y) && domain.is_periodic_bound(Y)) {
                //     jm = imod(j - 1, ny);
                //     jp = imod(j + 1, ny);
                // } else if (dim == static_cast<int>(Axis::Z) &&
                //            domain.is_periodic_bound(Z)) {
                //     km = imod(k - 1, nz);
                //     kp = imod(k + 1, nz);
                // }

                for (int ax = 0; ax < COMPONENTS; ++ax) {
                    const int row = vind(i, j, k, ax);

                    int col_m, col_c, col_p;
                    if (dim == static_cast<int>(Axis::X)) {
                        col_m = vind(im, j, k, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(ip, j, k, ax);
                    } else if (dim == static_cast<int>(Axis::Y)) {
                        col_m = vind(i, jm, k, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(i, jp, k, ax);
                    } else {
                        col_m = vind(i, j, km, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(i, j, kp, ax);
                    }

                    addTriplet(trips, domain, row, col_m, 0.25);
                    addTriplet(trips, domain, row, col_c, 0.5);
                    addTriplet(trips, domain, row, col_p, 0.25);
                }
            }
        }
    }

    mat.setFromTriplets(trips.begin(), trips.end());
}
