#ifndef MESH_H_
#define MESH_H_
#include "World.h"
#include "Read.h"
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include "bmatrix.h"

double3 interpolateE_linear(const Field3d& fieldE,
                            const double3& normalized_coord);
double3 interpolateE_Chen(const Field3d& fieldE, const double3& coord,
                          const Domain& domain);
double3 interpolateE(const Field3d& fieldE, const double3& normalized_coord,
                     ShapeType type);
double3 get_fieldB_in_pos(const Field3d& fieldB, const double3& coord,
                          const Domain& domain);
double3 get_fieldB_in_pos_new(const Field3d& field, const double3& coord,
                              const Domain& domain);
double3 interpolateE_ngp(const Field3d& fieldE,
                         const double3& normalized_coord);
double3 interpolateB_ngp(const Field3d& fieldB,
                         const double3& normalized_coord);
double calc_JE(const Field3d& fieldE, const Field3d& fieldJ,
               const Bounds& bounds);
double3 calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ,
                          const Bounds& bounds);
void apply_periodic_border_with_add(Field3d& field, const Bounds& bounds);

// void get_fields_in_pos(const Field3d& fieldE,const Field3d& fieldB, const
// double3& r, double3& locE, double3 &locB);

struct IndexingParams {
    int size_i, size_j, size_k;
    std::function<int(int, int, int)> index_func;
    int direction;
    struct Offsets {
        int i, j, k;
    } row_offs, col_offs;
};


struct Mesh{
    Mesh(){};
    void init(const Domain& domain,
              const ParametersMap& parameters);

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

    void stencil_curlE_openZ(Operator& mat, const Domain& domain);
    void stencil_curlB_openZ(Operator& mat, const Domain& domain);
    void apply_open_boundaries_z(std::vector<IndexMap>& LmatX);
    void apply_open_boundaries_z(Field3d& field);
    void print_operator(const Operator& oper);
    void processDirection(std::vector<Trip>& trips, const Block& block,
                                const IndexingParams& row_params,
                                const IndexingParams& col_params,
                                int output_direction, int i_cell, int j_cell,
                                int k_cell, double tolerance);
    std::vector<IndexMap> LmatX;

    //Sources and fields on the grid
    Field3d chargeDensityOld;
    Field3d chargeDensity;
    Operator divE;
    Bounds bounds;
    // Field multiply_LmatX2_vector(BMatrix2& LmatX2, const Field& vec,
    //                              const Domain& domain);

    inline int sind(int i, int j, int k) const {
        return i * ySize * zSize + j * zSize + k;
    };
    // index for 3D vector fields
    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return d + nd*(i * ySize * zSize + j * zSize + k);
    };
    // void convert_LmatX2_to_CSR(BMatrix2& LmatX2,
    //                                   const Domain& domain,
    //                                   MatrixCSR& csr_matrix);
    inline int pos_vind(int index, int n){
        std::vector<int> dim = {xSize, ySize, zSize, 3};
        int capacity = 1;
        for(unsigned int i = n + 1; i < dim.size(); i++){
            capacity *= dim[i];
        }
        return (index / capacity) % dim[n];
    }
    void zeroBoundL(Operator & mat) ;
    void zeroBoundJ(Field3d& field);
    void prepare();
    void computeB(const Field3d& fieldE, const Field3d& fieldEn,
                  Field3d& fieldB, double dt);
    void fdtd_explicit(Field3d& E, Field3d& B,
                       const Field3d& J, const double dt);
    void set_mirrors();
    void update_Lmat(const double3& coord, const Domain& domain, double charge,
                     double mass, double mpw, const Field3d& fieldB, const double dt);
    void update_LmatNGP(const double3& coord, const Domain& domain, double charge,
                      double mass, double mpw, const Field3d& fieldB,
                      const double dt);

    void update_Lmat2(const double3& coord, const Domain& domain, double charge,
                     double mass, double mpw, const Field3d& fieldB,
                     const double dt);

    void update_Lmat2_NGP(const double3& coord, const Domain& domain, double charge,
                     double mass, double mpw, const Field3d& fieldB,
                     const double dt);

    void apply_periodic_boundaries(std::vector<IndexMap>& LmatX);
    void apply_open_boundaries(std::vector<IndexMap>& LmatX,
                                 const Domain& domain);
    void apply_boundaries(std::vector<IndexMap>& LmatX, const Domain& domain);

    void apply_periodic_boundaries(Operator& LmatX);
    void apply_open_boundaries(Operator& LmatX, Domain& domain);
    void apply_boundaries(Operator& LmatX, Domain& domain);
    
    void apply_periodic_boundaries(Field3d& field);
    void apply_open_boundaries(Field3d& field, const Domain& domain);
    void apply_boundaries(Field3d& field, const Domain& domain);

    void apply_density_periodic_boundaries(Field3d& field);
    void apply_density_open_boundaries(Field3d& field, const Domain& domain);
    void apply_density_boundaries(Field3d& field, const Domain& domain);

    void set_uniform_field(Field3d& field, double bx, double by, double bz);

    double3 get_fieldE_in_cell(const Field3d& fieldE, int i, int j,
                               int k) const;
    double3 get_fieldB_in_cell(const Field3d& fieldB, int i, int j,
                               int k) const;
    double calc_energy_field(const Field3d& field) const;

    // double get_fieldE_energy() const { return calc_energy_field(fieldE); };
    // double get_fieldB_energy() const{
    //     return calc_energy_field(fieldB); //-calc_energy_field(fieldB0);
    // };

    ~Mesh(){
    }

    void stencil_curlB(Operator& mat, const Domain& domain);
    void stencil_curlE(Operator& mat, const Domain& domain);
    void stencil_Imat(Operator& mat, const Domain& domain);

    void stencil_Lmat(Operator& mat, const Domain& domain);
    void stencil_Lmat2(Operator& mat, const Domain& domain);
    void stencil_Lmat2_NGP(Operator& mat, const Domain& domain);
    template <typename IndexerX, typename IndexerY, typename IndexerZ,
          typename MatrixType>
    void convert_block_to_crs_format(MatrixType bmatrix, Operator& mat,
                                       const Domain& domain);
    void stencil_divE(Operator& mat, const Domain& domain);
    void stencil_curlE_periodic(std::vector<Trip>& trips, const Domain& domain);
   // void stencil_curlE_openZ(std::vector<Trip>& trips, const Domain& domain);
    void stencil_curlB_periodic(std::vector<Trip>& trips, const Domain& domain);
   // void stencil_curlB_openZ(std::vector<Trip>& trips, const Domain& domain);

    void predictE(Field3d& Ep, const Field3d& E, const Field3d& B,
                  Field3d& J, const double dt);
    void correctE(Field3d& En, const Field3d& E, const Field3d& B,
                  Field3d& J, const double dt);
    void impicit_find_fieldE(Field3d& Enew, const Field3d& E, const Field3d& B,
                             const Field3d& J, const double dt);
    double calculate_residual(const Field3d& Enew, const Field3d& E, const Field3d& B,
                             const Field3d& J, const double dt);
    void compute_fieldB(Field3d& Bn, const Field3d& B, const Field3d& E,
                        const Field3d& En, double dt);

    // general indexing routine (row major)
    inline constexpr int ind2(int x, int y, int z, int Nx, int Ny, int Nz) {
        return z + Nz * (y + Ny * x) + Nx*0;
    }


    void print_Lmat(const Operator& oper) {
        for (int i = 0; i < oper.outerSize(); ++i) {
            auto ix = pos_vind(i, 0);
            auto iy = pos_vind(i, 1);
            auto iz = pos_vind(i, 2);
            auto id = pos_vind(i, 3);
            for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(oper, i);
                 it; ++it) {
                auto ind2 = it.col();
                auto value = it.value();
                auto ix1 = pos_vind(ind2, 0);
                auto iy1 = pos_vind(ind2, 1);
                auto iz1 = pos_vind(ind2, 2);
                auto id1 = pos_vind(ind2, 3);
                //if (ind2 < 0 || ind2 >= oper.outerSize()) {
                if (value != 0)
                    std::cout << ix << " " << iy << " " << iz << " " << id
                              << " " << ix1 << " " << iy1 << " " << iz1 << " "
                              << id1 << " " << value << "\n";
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

void print_operator(const Operator &oper);

#include "operators.h"

#endif 
