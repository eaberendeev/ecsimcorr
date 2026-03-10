// #include "World.h"
// #include "operators.h"
// #include "util.h"
// // matrix RowMajor
// // (0,0) (0,1) (0,2)
// // (1,0) (1,1) (1,2)
// // (2,0) (2,1) (2,2)

// // Ex(i+/2,j,k), Ey(i,j+1/2,k), Ez(i,j,k+1/2)
// // Bx(i,j+1/2,k+1/2), By(i+1/2,j,k+1/2), Bz(i+/2,j+1/2,k)

// inline void stencil_Imat(Operator &mat, const Domain &domain) {
//     std::vector<Trip> trips;
//     const auto size = domain.size();
//     int totalSize = size.x() * size.y() * size.z();
//     trips.reserve(totalSize);
//     auto vind =
//         [&](int i, int j, int k, int d) { return domain.vind(i, j, k, d); };

//     for (int i = 0; i < size.x(); i++) {
//         for (int j = 0; j < size.y(); j++) {
//             for (int k = 0; k < size.z(); k++) {
//                 trips.push_back(Trip(vind(i, j, k, 0), vind(i, j, k, 0), 1.0));
//                 trips.push_back(Trip(vind(i, j, k, 1), vind(i, j, k, 1), 1.0));
//                 trips.push_back(Trip(vind(i, j, k, 2), vind(i, j, k, 2), 1.0));
//             }
//         }
//     }

//     mat.setFromTriplets(trips.begin(), trips.end());
//     std::cout << "Identity operator has been created succsessfully.\n";
//     std::cout << "rows: " << mat.rows() << " cols: " << mat.cols()
//               << " elements: " << trips.size() << std::endl;
// }

// inline void stencil_curlB(Operator &mat, const Domain &domain) {
//     std::vector<Trip> trips;
//     const auto size = domain.size();
//     int totalSize = size.x() * size.y() * size.z() * 12;
//     trips.reserve(totalSize);
//     const double dx = domain.cell_size().x();
//     const double dy = domain.cell_size().y();
//     const double dz = domain.cell_size().z();
//     auto vind =
//         [&](int i, int j, int k, int d) { return domain.vind(i, j, k, d); };

//     auto addTriplet = [](std::vector<Trip> &trips, const Domain &domain,
//                          int vindElec, int vindMag, double val) {
//         bool onArea = domain.in_region_electric(vindElec);
//         if (onArea) {
//             trips.push_back(Trip(vindElec, vindMag, val));
//         }
//     };

//     for (int i = 0; i < size.x(); i++) {
//         for (int j = 0; j < size.y(); j++) {
//             for (int k = 0; k < size.z(); k++) {
//                 const int im = (i != 0) || !domain.is_periodic_bound(X)
//                                    ? i - 1
//                                    : size.x() - 4;
//                 const int jm = (j != 0) || !domain.is_periodic_bound(Y)
//                                    ? j - 1
//                                    : size.y() - 4;
//                 const int km = (k != 0) || !domain.is_periodic_bound(Z)
//                                    ? k - 1
//                                    : size.z() - 4;

//                 const int vindx = vind(i, j, k, 0);
//                 const int vindy = vind(i, j, k, 1);
//                 const int vindz = vind(i, j, k, 2);

//                 // (x)[i+1/2,j,k]
//                 // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
//                 double val = 1.0 / dy;
//                 addTriplet(trips, domain, vindx, vind(i, j, k, 2), val);
//                 addTriplet(trips, domain, vindx, vind(i, jm, k, 2), -val);
//                 // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
//                 val = -1.0 / dz;
//                 addTriplet(trips, domain, vindx, vind(i, j, k, 1), val);
//                 addTriplet(trips, domain, vindx, vind(i, j, km, 1), -val);

//                 // (y)[i,j+1/2,k]
//                 // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
//                 val = 1.0 / dz;
//                 addTriplet(trips, domain, vindy, vind(i, j, k, 0), val);
//                 addTriplet(trips, domain, vindy, vind(i, j, km, 0), -val);
//                 // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
//                 val = -1.0 / dx;
//                 addTriplet(trips, domain, vindy, vind(i, j, k, 2), val);
//                 addTriplet(trips, domain, vindy, vind(im, j, k, 2), -val);

//                 // (z)[i,j,k+1/2]
//                 // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
//                 val = 1.0 / dx;
//                 addTriplet(trips, domain, vindz, vind(i, j, k, 1), val);
//                 addTriplet(trips, domain, vindz, vind(im, j, k, 1), -val);
//                 // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
//                 val = -1.0 / dy;
//                 addTriplet(trips, domain, vindz, vind(i, j, k, 0), val);
//                 addTriplet(trips, domain, vindz, vind(i, jm, k, 0), -val);
//             }
//         }
//     }
//     mat.setFromTriplets(trips.begin(), trips.end());
//     std::cout << "CurlE operator has been created succsessfully.\n";
//     std::cout << "rows: " << mat.rows() << " cols: " << mat.cols()
//               << " elements: " << trips.size() << std::endl;
// }

// inline void stencil_curlE(Operator &mat, const Domain &domain) {
//     std::vector<Trip> trips;
//     const auto size = domain.size();
//     int totalSize = size.x() * size.y() * size.z() * 12;
//     trips.reserve(totalSize);
//     double dx = domain.cell_size().x();
//     double dy = domain.cell_size().y();
//     double dz = domain.cell_size().z();

//     auto vind =
//         [&](int i, int j, int k, int d) { return domain.vind(i, j, k, d); };

//     auto addTriplet = [](std::vector<Trip> &trips, const Domain &domain,
//                          int vindMag, int vindElec, double val) {
//         bool onArea = domain.in_region_magnetic(vindMag);
//         if (onArea) {
//             trips.push_back(Trip(vindMag, vindElec, val));
//         }
//     };

//     for (int i = 0; i < size.x(); i++) {
//         for (int j = 0; j < size.y(); j++) {
//             for (int k = 0; k < size.z(); k++) {
//                 const int ip =
//                     (i != size.x() - 1) || !domain.is_periodic_bound(X) ? i + 1
//                                                                         : 3;
//                 const int jp =
//                     (j != size.y() - 1) || !domain.is_periodic_bound(Y) ? j + 1
//                                                                         : 3;
//                 const int kp =
//                     (k != size.z() - 1) || !domain.is_periodic_bound(Z) ? k + 1
//                                                                         : 3;
//                 const int vindx = vind(i, j, k, 0);
//                 const int vindy = vind(i, j, k, 1);
//                 const int vindz = vind(i, j, k, 2);

//                 // (x)[i,j+1/2,k+1/2]
//                 // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
//                 double val = 1.0 / dy;
//                 addTriplet(trips, domain, vindx, vind(i, jp, k, 2), val);
//                 addTriplet(trips, domain, vindx, vind(i, j, k, 2), -val);
//                 // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
//                 val = -1.0 / dz;
//                 addTriplet(trips, domain, vindx, vind(i, j, kp, 1), val);
//                 addTriplet(trips, domain, vindx, vind(i, j, k, 1), -val);

//                 // (y)[i+1/2,j,k+1/2]
//                 // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
//                 val = 1.0 / dz;
//                 addTriplet(trips, domain, vindy, vind(i, j, kp, 0), val);
//                 addTriplet(trips, domain, vindy, vind(i, j, k, 0), -val);
//                 // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
//                 val = -1.0 / dx;
//                 addTriplet(trips, domain, vindy, vind(ip, j, k, 2), val);
//                 addTriplet(trips, domain, vindy, vind(i, j, k, 2), -val);

//                 // (z)[i+1/2,j+1/2,k]
//                 // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
//                 val = 1.0 / dx;
//                 addTriplet(trips, domain, vindz, vind(ip, j, k, 1), val);
//                 addTriplet(trips, domain, vindz, vind(i, j, k, 1), -val);
//                 // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
//                 val = -1.0 / dy;
//                 addTriplet(trips, domain, vindz, vind(i, jp, k, 0), val);
//                 addTriplet(trips, domain, vindz, vind(i, j, k, 0), -val);
//             }
//         }
//     }
//     mat.setFromTriplets(trips.begin(), trips.end());
//     std::cout << "CurlB operator has been created succsessfully.\n";
//     std::cout << "rows: " << mat.rows() << " cols: " << mat.cols()
//               << " elements: " << trips.size() << std::endl;
// }

// inline void stencil_divE(Operator &mat, const Domain &domain) {
//     // !!!!! needs bound condition and if cases!!!!!!
//     std::vector<Trip> trips;
//     const auto size = domain.size();
//     int totalSize = size.x() * size.y() * size.z() * 6;
//     trips.reserve(totalSize);
//     double dx = domain.cell_size().x();
//     double dy = domain.cell_size().y();
//     double dz = domain.cell_size().z();
//     auto vind =
//         [&](int i, int j, int k, int d) { return domain.vind(i, j, k, d); };

//     for (int i = 0; i < size.x(); i++) {
//         for (int j = 0; j < size.y(); j++) {
//             for (int k = 0; k < size.z(); k++) {
//                 const int im = (i != 0) ? i - 1 : size.x() - 4;
//                 const int jm = (j != 0) ? j - 1 : size.y() - 4;
//                 const int km = (k != 0) ? k - 1 : size.z() - 4;

//                 const int sindx = domain.sind(i, j, k);

//                 // [i,j,k]
//                 // ( Ex[i+1/2,j,k] - Ex[i-1,j,k] ) / dx
//                 double val = 1.0 / dx;
//                 trips.push_back(Trip(sindx, vind(i, j, k, 0), val));
//                 trips.push_back(Trip(sindx, vind(im, j, k, 0), -val));
//                 // ( Ex[i,j+1/2,k] - Ex[i,j-1/2,k] ) / dy
//                 val = 1.0 / dy;
//                 trips.push_back(Trip(sindx, vind(i, j, k, 1), val));
//                 trips.push_back(Trip(sindx, vind(i, jm, k, 1), -val));
//                 // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
//                 val = 1.0 / dz;
//                 trips.push_back(Trip(sindx, vind(i, j, km, 2), val));
//             }
//         }
//     }
//     mat.setFromTriplets(trips.begin(), trips.end());
//     std::cout << "DivE operator has been created succsessfully.\n";
//     std::cout << "rows: " << mat.rows() << " cols: " << mat.cols()
//               << " elements: " << trips.size() << std::endl;
// }
