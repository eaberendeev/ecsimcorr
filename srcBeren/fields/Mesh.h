#ifndef MESH_H_
#define MESH_H_
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "Read.h"
#include "World.h"
#include "bmatrix.h"

struct Mesh {
    Mesh(){};
    void init(const Domain& domain, double dt);

    void stencil_smooth_1d(Operator& mat, const Domain& domain, int dim);
    Operator Lmat;
    Operator Lmat2;
    Operator Mmat;
    Operator Imat;
    Operator curlE;
    Operator curlB;
    Operator IMmat;

    BlockMatrix LmatX2;
    BlockMatrixNGP LmatX_NGP;

    void print_operator(const Operator& oper);

    // Sources and fields on the grid
    Field3d chargeDensityOld;
    Field3d chargeDensity;
    Operator divE;

    inline int sind(int i, int j, int k) const {
        return i * ySize * zSize + j * zSize + k;
    };
    // index for 3D vector fields
    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return d + nd * (i * ySize * zSize + j * zSize + k);
    };

    inline int pos_vind(int index, int n) {
        std::vector<int> dim = {xSize, ySize, zSize, 3};
        int capacity = 1;
        for (unsigned int i = n + 1; i < dim.size(); i++) {
            capacity *= dim[i];
        }
        return (index / capacity) % dim[n];
    }

    void prepare();
    void computeB(const Field3d& fieldE, const Field3d& fieldEn, Field3d& fieldB, double dt);
    void fdtd_explicit(Field3d& E, Field3d& B, const Field3d& J, const double dt);
    void set_mirrors();
    void update_Lmat(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                     const Field3d& fieldB, const double dt);
    void update_LmatNGP(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                        const Field3d& fieldB, const double dt);

    void update_Lmat2(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                      const Field3d& fieldB, const double dt);

    void update_Lmat2_NGP(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                          const Field3d& fieldB, const double dt);

    void set_uniform_field(Field3d& field, double bx, double by, double bz);

    double calc_energy_field(const Field3d& field) const;

    ~Mesh() {
    }

    void stencil_curlB(Operator& mat, const Domain& domain);
    void stencil_curlE(Operator& mat, const Domain& domain);
    void stencil_Imat(Operator& mat, const Domain& domain);

    void stencil_Lmat(Operator& mat, const Domain& domain);
    void stencil_Lmat2(Operator& mat, const Domain& domain);
    void stencil_Lmat2_NGP(Operator& mat, const Domain& domain);
    template <typename IndexerX, typename IndexerY, typename IndexerZ, typename MatrixType>
    void convert_block_to_crs_format(MatrixType bmatrix, Operator& mat, const Domain& domain);
    void stencil_divE(Operator& mat, const Domain& domain);
    void stencil_curlE_periodic(std::vector<Trip>& trips, const Domain& domain);
    void stencil_curlB_periodic(std::vector<Trip>& trips, const Domain& domain);

    void impicit_find_fieldE(Field3d& Enew, const Field3d& E, const Field3d& B, const Field3d& J, const double dt);
    double calculate_residual(const Field3d& Enew, const Field3d& E, const Field3d& B, const Field3d& J,
                              const double dt);
    void compute_fieldB(Field3d& Bn, const Field3d& B, const Field3d& E, const Field3d& En, double dt);

    // general indexing routine (row major)
    inline constexpr int ind2(int x, int y, int z, int Nx, int Ny, int Nz) {
        return z + Nz * (y + Ny * x) + Nx * 0;
    }

    void print_Lmat(const Operator& oper) {
        for (int i = 0; i < oper.outerSize(); ++i) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);
            for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(oper, i); it; ++it) {
                auto ind2 = it.col();
                auto value = it.value();
                auto ix1 = pos_vind(ind2, 0);
                auto iy1 = pos_vind(ind2, 1);
                auto iz1 = pos_vind(ind2, 2);
                auto id1 = pos_vind(ind2, 3);
                // if (ind2 < 0 || ind2 >= oper.outerSize()) {
                if (value != 0)
                    std::cout << ix << " " << iy << " " << iz << " " << id << " " << ix1 << " " << iy1 << " " << iz1
                              << " " << id1 << " " << value << "\n";
                //}
            }
        }
    }

   private:
    double xCellSize;
    double yCellSize;
    double zCellSize;
    int xSize;
    int ySize;
    int zSize;
};

void print_operator(const Operator& oper);

#include "operators.h"

#endif
