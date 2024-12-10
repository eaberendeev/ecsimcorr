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

 void Mesh::stencil_Imat(BMatrix& mat, const Domain &domain)
 { 
 // !!!!! needs bound condition and if cases!!!!!!
  const auto size = domain.size();
    
  for(int i = 0; i < size.x(); i++){
    for(int j = 0; j < size.y(); j++){
      for(int k = 0; k < size.z(); k++){
        // i,j,k
        mat[vind(i, j, k, X)][vind(i, j, k, X)] = 1.0;
        mat[vind(i, j, k, Y)][vind(i, j, k, Y)] = 1.0;
        mat[vind(i, j, k, Z)][vind(i, j, k, Z)] = 1.0;
      }
    }
  }
}

void Mesh::stencil_curlB(BMatrix &mat, const Domain &domain) {
    // TO DO: create a different boundary cases
    // NOW X and Y always periodic
    if (domain.is_periodic_bound(Z)) {
        stencil_curlB_periodic(mat, domain);
    } else {
        stencil_curlB_openZ(mat, domain);
    }
}

void Mesh::stencil_curlB_periodic(BMatrix &mat, const Domain &domain) {
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
                mat[vindx][vind(i, j, k, Z)] = val;
                mat[vindx][vind(i, jm, k, Z)] = -val;
                // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                val = -1.0 / dz;
                mat[vindx][vind(i, j, k, Y)] = val;
                mat[vindx][vind(i, j, km, Y)] = -val;

                // (y)[i,j+1/2,k]
                // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                val = 1.0 / dz;
                mat[vindy][vind(i, j, k, X)] = val;
                mat[vindy][vind(i, j, km, X)] = -val;

                // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                val = -1.0 / dx;
                mat[vindy][vind(i, j, k, Z)] = val;
                mat[vindy][vind(im, j, k, Z)] = -val;
                // (z)[i,j,k+1/2]
                // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                val = 1.0 / dx;
                mat[vindz][vind(i, j, k, Y)] = val;
                mat[vindz][vind(im, j, k, Y)] = -val;
                // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                val = -1.0 / dy;
                mat[vindz][vind(i, j, k, X)] = val;
                mat[vindz][vind(i, jm, k, X)] = -val;
            }
        }
    }
}

// TODO merge with periodic
void Mesh::stencil_curlB_openZ(BMatrix &mat, const Domain &domain) {
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
                        mat[vindz][vind(i, j, k, Y)] = val;
                        mat[vindz][vind(im, j, k, Y)] = -val;

                        // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                        val = -1.0 / dy;
                        mat[vindz][vind(i, j, k, X)] = val;
                        mat[vindz][vind(i, jm, k, X)] = -val;
                    }
                    continue;
                }

                // (x)[i+1/2,j,k]
                // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                double val = 1.0 / dy;
                mat[vindx][vind(i, j, k, Z)] = val;
                mat[vindx][vind(i, jm, k, Z)] = -val;
                // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                val = -1.0 / dz;
                mat[vindx][vind(i, j, k, Y)] = val;
                mat[vindx][vind(i, j, km, Y)] = -val;
                // (y)[i,j+1/2,k]
                // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                val = 1.0 / dz;
                mat[vindy][vind(i, j, k, X)] = val;
                mat[vindy][vind(i, j, km, X)] = -val;
                // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                val = -1.0 / dx;
                mat[vindy][vind(i, j, k, Z)] = val;
                mat[vindy][vind(im, j, k, Z)] = -val;
                // (z)[i,j,k+1/2]
                // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                val = 1.0 / dx;
                mat[vindz][vind(i, j, k, Y)] = val;
                mat[vindz][vind(im, j, k, Y)] = -val;
                // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                val = -1.0 / dy;
                mat[vindz][vind(i, j, k, X)] = val;
                mat[vindz][vind(i, jm, k, X)] = -val;
            }
        }
    }
}

void Mesh::stencil_curlE(BMatrix &mat, const Domain &domain) {
    // TO DO: create a different boundary cases
    // NOW X and Y always periodic
    if (domain.is_periodic_bound(Z)) {
        stencil_curlE_periodic(mat, domain);
    } else {
        stencil_curlE_openZ(mat, domain);
    }
}

void Mesh::stencil_curlE_periodic(BMatrix &mat, const Domain &domain) {
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
                mat[vindx][vind(i, jp, k, Z)] = val;
                mat[vindx][vind(i, j, k, Z)] = -val;
                // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                val = -1.0 / dz;
                mat[vindx][vind(i, j, kp, Y)] = val;
                mat[vindx][vind(i, j, k, Y)] = -val;
                // (y)[i+1/2,j,k+1/2]
                // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                val = 1.0 / dz;
                mat[vindy][vind(i, j, kp, X)] = val;
                mat[vindy][vind(i, j, k, X)] = -val;
                // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                val = -1.0 / dx;
                mat[vindy][vind(ip, j, k, Z)] = val;
                mat[vindy][vind(i, j, k, Z)] = -val;

                // (z)[i+1/2,j+1/2,k]
                // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                val = 1.0 / dx;
                mat[vindz][vind(ip, j, k, Y)] = val;
                mat[vindz][vind(i, j, k, Y)] = -val;
                // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                val = -1.0 / dy;
                mat[vindz][vind(i, jp, k, X)] = val;
                mat[vindz][vind(i, j, k, X)] = -val;
            }
        }
    }
}

void Mesh::stencil_curlE_openZ(BMatrix &mat, const Domain &domain) {
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
                    mat[vindx][vind(i, jp, k, Z)] = val;
                    mat[vindx][vind(i, j, k, Z)] = -val;
                }
                // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                val = -1.0 / dz;
                if (kp > 1 && kp < size.z() - 2) {
                    mat[vindx][vind(i, j, kp, Y)] = val;
                }
                if (k > 1 && k < size.z() - 2) {
                    mat[vindx][vind(i, j, k, Y)] = -val;
                }
                // (y)[i+1/2,j,k+1/2]
                // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                val = 1.0 / dz;
                if (kp > 1 && kp < size.z() - 2) {
                    mat[vindy][vind(i, j, kp, X)] = val;
                }
                if (k > 1 && k < size.z() - 2) {
                    mat[vindy][vind(i, j, k, X)] = -val;
                }
                // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                val = -1.0 / dx;
                if (k > 0 && k < size.z() - 2) {
                    mat[vindy][vind(ip, j, k, Z)] = val;
                    mat[vindy][vind(i, j, k, Z)] = -val;
                }
                // (z)[i+1/2,j+1/2,k]
                // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                if (k > 1 && k < size.z() - 2) {
                    val = 1.0 / dx;
                    mat[vindz][vind(ip, j, k, Y)] = val;
                    mat[vindz][vind(i, j, k, Y)] = -val;
                    // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                    val = -1.0 / dy;
                    mat[vindz][vind(i, jp, k, X)] = val;
                    mat[vindz][vind(i, j, k, X)] = -val;
                }
            }
        }
    }
}

void Mesh::stencil_divE(BMatrix &mat, const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    const auto size = domain.size();
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();
    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int im = (i != 0) ? i - 1 : size.x() - 4;
                const int jm = (j != 0) ? j - 1 : size.y() - 4;
                const int km = (k != 0) ? k - 1 : size.z() - 4;

                const int sindx = sind(i, j, k);

                // [i,j,k]
                // ( Ex[i+1/2,j,k] - Ex[i-1,j,k] ) / dx
                double val = 1.0 / dx;
                mat[sindx][vind(i, j, k, X)] = val;
                mat[sindx][vind(im, j, k, X)] = -val;
                // ( Ex[i,j+1/2,k] - Ex[i,j-1/2,k] ) / dy
                val = 1.0 / dy;
                mat[sindx][vind(i, j, k, Y)] = val;
                mat[sindx][vind(i, jm, k, Y)] = -val;
                // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
                val = 1.0 / dz;
                mat[sindx][vind(i, j, k, Z)] = val;
                mat[sindx][vind(i, j, km, Z)] = -val;
            }
        }
    }
}
