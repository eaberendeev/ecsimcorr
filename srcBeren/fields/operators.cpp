#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "util.h"

// matrix ColMajor
// (0,0) (0,1) (0,2)
// (1,0) (1,1) (1,2)
// (2,0) (2,1) (2,2)

// Ex(i+/2,j,k), Ey(i,j+1/2,k), Ez(i,j,k+1/2) 
// Bx(i,j+1/2,k+1/2), By(i+1/2,j,k+1/2), Bz(i+/2,j+1/2,k)

 void Mesh::stencil_Imat(const Domain &domain)
 { 
 // !!!!! needs bound condition and if cases!!!!!!
  std::vector<Trip> trips;
  const auto size = domain.size();
  int totalSize = size.x()*size.y()*size.z();
  trips.reserve(totalSize);
    
  for(int i = 0; i < size.x(); i++){
    for(int j = 0; j < size.y(); j++){
      for(int k = 0; k < size.z(); k++){

        // i,j,k
        trips.push_back(Trip(vind(i,j,k,0),vind(i,j,k,0),1.0));

        // i,j,k
        trips.push_back(Trip(vind(i,j,k,1),vind(i,j,k,1),1.0));
          
        // i,j,k
        trips.push_back(Trip(vind(i,j,k,2),vind(i,j,k,2),1.0));
          
      }
    }
  }
  std::cout << size << " " << 3*totalSize << " " << Imat.rows() << " " << Imat.cols() <<  "\n";
  Imat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_Lmat(const Domain &domain) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    const int rowsCount = 3 * size.x() * size.y() * size.z();
    const size_t totalSize = (size_t)rowsCount * LMAT_MAX_ELEMENTS_PER_ROW;
    std::cout << totalSize << "\n";
    trips.reserve(totalSize);

    for (int row = 0; row < rowsCount; row++) {
        for (const auto &[col, value] : LmatX[row]) {
            if (std::abs(value) > LMAT_VALUE_TOLERANCE)
                trips.emplace_back(row, col, value);
        }
    }

    Lmat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_Lmat2(const Domain &domain) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    const int rowsCount = 3 * size.x() * size.y() * size.z();
    const size_t totalSize = (size_t) rowsCount * LMAT_MAX_ELEMENTS_PER_ROW;
    std::cout << totalSize << "\n";
    trips.reserve(totalSize);

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                int row = vind(i, j, k, 0);
                for (int i1 = 0; i1 < 2; i1++) {
                    for (int j1 = 0; j1 < 2; j1++) {
                        for (int k1 = 0; k1 < 2; k1++) {
                            int col = vind(i + i1, j + j1, k + k1, 0);
                            double value = LmatX2[row][get_col_index_Lx(i1, j1, k1, X)];
                            if (std::abs(value) > LMAT_VALUE_TOLERANCE){
                                std::cout << get_col_index_Lx(i1, j1, k1, X)
                                           << " " << value << "\n";
                                    trips.emplace_back(row, col, value);
                            }
                        }
                    }
                }
            }
        }
    }
    Lmat2.setFromTriplets(trips.begin(), trips.end());
}
// TO DO: check if this is correct

// void Mesh::stencil_Lmat(const Domain &domain) {
//     const auto size = domain.size();
//     const int totalSize =
//         3 * size.x() * size.y() * size.z() * LMAT_MAX_ELEMENTS_PER_ROW;
//     const double tolerance = 1.e-16;

//     std::vector<std::vector<Trip>> thread_trips(omp_get_max_threads());

// #pragma omp parallel
//     {
//         int thread_num = omp_get_thread_num();
//         std::vector<Trip> &trips = thread_trips[thread_num];
//         trips.reserve(totalSize / omp_get_max_threads());

// #pragma omp for
//         for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
//             for (auto it = LmatX[i].begin(); it != LmatX[i].end(); ++it) {
//                 int ind2 = it->first;
//                 double value = it->second;
//                 if (std::abs(value) > tolerance) {
//                     trips.emplace_back(i, ind2, value);
//                 }
//             }
//         }
//     }

//     std::vector<Trip> merged_trips;
//     merged_trips.reserve(totalSize);

//     for (const auto &thread_trip : thread_trips) {
//         merged_trips.insert(merged_trips.end(), thread_trip.begin(),
//                             thread_trip.end());
//     }

//     Lmat.setFromTriplets(merged_trips.begin(), merged_trips.end());
// }

void Mesh::stencil_curlB(const Domain &domain) {
  // TO DO: create a different boundary cases
  // NOW X and Y always periodic
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x()*size.y()*size.z()*12;
    trips.reserve(totalSize);
    if (domain.is_periodic_bound(Z)) {
        stencil_curlB_periodic(trips, domain);
    } else {
        stencil_curlB_openZ(trips, domain);
    }

    curlB.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlB_periodic(std::vector<Trip> &trips,
                                  const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    const auto size = domain.size();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int im = (i != 0) ? i - 1 : size.x() - 4;
                const int jm = (j != 0) ? j - 1 : size.y() - 4;
                const int km = (k != 0) ? k - 1 : size.z() - 4;

                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                // (x)[i+1/2,j,k]
                // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                double val = 1.0 / dy;
                trips.push_back(Trip(vindx, vind(i, j, k, 2), val));
                trips.push_back(Trip(vindx, vind(i, jm, k, 2), -val));
                // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                val = -1.0 / dz;
                trips.push_back(Trip(vindx, vind(i, j, k, 1), val));
                trips.push_back(Trip(vindx, vind(i, j, km, 1), -val));

                // (y)[i,j+1/2,k]
                // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                val = 1.0 / dz;
                trips.push_back(Trip(vindy, vind(i, j, k, 0), val));
                trips.push_back(Trip(vindy, vind(i, j, km, 0), -val));
                // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                val = -1.0 / dx;
                trips.push_back(Trip(vindy, vind(i, j, k, 2), val));
                trips.push_back(Trip(vindy, vind(im, j, k, 2), -val));

                // (z)[i,j,k+1/2]
                // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                val = 1.0 / dx;
                trips.push_back(Trip(vindz, vind(i, j, k, 1), val));
                trips.push_back(Trip(vindz, vind(im, j, k, 1), -val));
                // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                val = -1.0 / dy;
                trips.push_back(Trip(vindz, vind(i, j, k, 0), val));
                trips.push_back(Trip(vindz, vind(i, jm, k, 0), -val));
            }
        }
    }
}

void Mesh::stencil_curlB_openZ(std::vector<Trip> &trips, const Domain &domain) {
    const auto size = domain.size();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int im = (i != 0) ? i - 1 : size.x() - 4;
                const int jm = (j != 0) ? j - 1 : size.y() - 4;
                const int km = (k != 0) ? k - 1 : size.z() - 4;
                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                if (k < 2 || k > size.z() - 3) {
                    if (k == 1) {
                        // (z)[i,j,k+1/2]
                        // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                        double val = 1.0 / dx;
                        trips.push_back(Trip(vindz, vind(i, j, k, 1), val));
                        trips.push_back(Trip(vindz, vind(im, j, k, 1), -val));
                        // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                        val = -1.0 / dy;
                        trips.push_back(Trip(vindz, vind(i, j, k, 0), val));
                        trips.push_back(Trip(vindz, vind(i, jm, k, 0), -val));
                    }
                    continue;
                }

                // (x)[i+1/2,j,k]
                // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                double val = 1.0 / dy;
                trips.push_back(Trip(vindx, vind(i, j, k, 2), val));
                trips.push_back(Trip(vindx, vind(i, jm, k, 2), -val));
                // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                val = -1.0 / dz;
                trips.push_back(Trip(vindx, vind(i, j, k, 1), val));
                trips.push_back(Trip(vindx, vind(i, j, km, 1), -val));

                // (y)[i,j+1/2,k]
                // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                val = 1.0 / dz;
                trips.push_back(Trip(vindy, vind(i, j, k, 0), val));
                trips.push_back(Trip(vindy, vind(i, j, km, 0), -val));
                // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                val = -1.0 / dx;
                trips.push_back(Trip(vindy, vind(i, j, k, 2), val));
                trips.push_back(Trip(vindy, vind(im, j, k, 2), -val));

                // (z)[i,j,k+1/2]
                // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                val = 1.0 / dx;
                trips.push_back(Trip(vindz, vind(i, j, k, 1), val));
                trips.push_back(Trip(vindz, vind(im, j, k, 1), -val));
                // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                val = -1.0 / dy;
                trips.push_back(Trip(vindz, vind(i, j, k, 0), val));
                trips.push_back(Trip(vindz, vind(i, jm, k, 0), -val));
            }
        }
    }
}

void Mesh::stencil_curlE(const Domain &domain) {
    // TO DO: create a different boundary cases
    // NOW X and Y always periodic
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);
    if (domain.is_periodic_bound(Z)) {
        stencil_curlE_periodic(trips, domain);
    } else {
        stencil_curlE_openZ(trips, domain);
    }
    curlE.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlE_periodic(std::vector<Trip> &trips,
                                  const Domain &domain) {
    const auto size = domain.size();
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int ip = (i != size.x() - 1) ? i + 1 : 3;
                const int jp = (j != size.y() - 1) ? j + 1 : 3;
                const int kp = (k != size.z() - 1) ? k + 1 : 3;

                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                // (x)[i,j+1/2,k+1/2]
                // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
                double val = 1.0 / dy;
                trips.push_back(Trip(vindx, vind(i, jp, k, 2), val));
                trips.push_back(Trip(vindx, vind(i, j, k, 2), -val));
                // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                val = -1.0 / dz;
                trips.push_back(Trip(vindx, vind(i, j, kp, 1), val));
                trips.push_back(Trip(vindx, vind(i, j, k, 1), -val));

                // (y)[i+1/2,j,k+1/2]
                // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                val = 1.0 / dz;
                trips.push_back(Trip(vindy, vind(i, j, kp, 0), val));
                trips.push_back(Trip(vindy, vind(i, j, k, 0), -val));
                // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                val = -1.0 / dx;
                trips.push_back(Trip(vindy, vind(ip, j, k, 2), val));
                trips.push_back(Trip(vindy, vind(i, j, k, 2), -val));

                // (z)[i+1/2,j+1/2,k]
                // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                val = 1.0 / dx;
                trips.push_back(Trip(vindz, vind(ip, j, k, 1), val));
                trips.push_back(Trip(vindz, vind(i, j, k, 1), -val));
                // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                val = -1.0 / dy;
                trips.push_back(Trip(vindz, vind(i, jp, k, 0), val));
                trips.push_back(Trip(vindz, vind(i, j, k, 0), -val));
            }
        }
    }
}

void Mesh::stencil_curlE_openZ(std::vector<Trip> &trips, const Domain &domain) {
    const auto size = domain.size();

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z() - 2; k++) {
                const int ip = (i != size.x() - 1) ? i + 1 : 3;
                const int jp = (j != size.y() - 1) ? j + 1 : 3;
                const int kp = k + 1;

                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                // (x)[i,j+1/2,k+1/2]
                // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
                double val = 1.0 / dy;
                if (k > 0 && k < size.z() - 2) {
                    trips.push_back(Trip(vindx, vind(i, jp, k, 2), val));
                    trips.push_back(Trip(vindx, vind(i, j, k, 2), -val));
                }
                // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                val = -1.0 / dz;
                if (kp > 1 && kp < size.z() - 2) {
                    trips.push_back(Trip(vindx, vind(i, j, kp, 1), val));
                }
                if (k > 1 && k < size.z() - 2) {
                    trips.push_back(Trip(vindx, vind(i, j, k, 1), -val));
                }
                // (y)[i+1/2,j,k+1/2]
                // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                val = 1.0 / dz;
                if (kp > 1 && kp < size.z() - 2) {
                    trips.push_back(Trip(vindy, vind(i, j, kp, 0), val));
                }
                if (k > 1 && k < size.z() - 2) {
                    trips.push_back(Trip(vindy, vind(i, j, k, 0), -val));
                }
                // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                val = -1.0 / dx;
                if (k > 0 && k < size.z() - 2) {
                    trips.push_back(Trip(vindy, vind(ip, j, k, 2), val));
                    trips.push_back(Trip(vindy, vind(i, j, k, 2), -val));
                }
                // (z)[i+1/2,j+1/2,k]
                // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                if (k > 1 && k < size.z() - 2) {
                    val = 1.0 / dx;
                    trips.push_back(Trip(vindz, vind(ip, j, k, 1), val));
                    trips.push_back(Trip(vindz, vind(i, j, k, 1), -val));
                    // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                    val = -1.0 / dy;
                    trips.push_back(Trip(vindz, vind(i, jp, k, 0), val));
                    trips.push_back(Trip(vindz, vind(i, j, k, 0), -val));
                }
            }
        }
    }
}

void Mesh::stencil_divE(const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x()*size.y()*size.z()*6;
    trips.reserve(totalSize);
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();
    for(int i = 0; i < size.x(); i++){
      for(int j = 0; j < size.y(); j++){
        for(int k = 0; k < size.z(); k++){
            const int im = (i != 0) ? i - 1 : size.x() - 4;
            const int jm = (j != 0) ? j - 1 : size.y() - 4;
            const int km = (k != 0) ? k - 1 : size.z() - 4;

            const int sindx = sind(i, j, k);

            // [i,j,k]
            // ( Ex[i+1/2,j,k] - Ex[i-1,j,k] ) / dx
            double val = 1.0 / dx;
            trips.push_back(Trip(sindx, vind(i, j, k, 0), val));
            trips.push_back(Trip(sindx, vind(im, j, k, 0), -val));
            // ( Ex[i,j+1/2,k] - Ex[i,j-1/2,k] ) / dy
            val = 1.0 / dy;
            trips.push_back(Trip(sindx, vind(i, j, k, 1), val));
            trips.push_back(Trip(sindx, vind(i, jm, k, 1), -val));
            // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
            val = 1.0 / dz;
            trips.push_back(Trip(sindx, vind(i, j, k, 2), val));
            trips.push_back(Trip(sindx, vind(i, j, km, 2), -val));

        }
      }
    }
    divE.setFromTriplets(trips.begin(), trips.end());
}
