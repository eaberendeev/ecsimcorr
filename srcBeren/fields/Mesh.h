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

double3 get_fieldE_in_pos(const Field3d& fieldE, const double3& coord,
                          const Domain& domain);
double3 interpolateE_Chen(const Field3d& fieldE, const double3& coord,
                          const Domain& domain);

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

using Array44 = std::array<double, 44>;

struct Mesh{
    Mesh(){};
    void init(const Domain& domain,
              const ParametersMap& parameters);

    Operator Lmat;
    Operator Lmat2;
    Operator Mmat;
    Operator Imat;
    Operator curlE;
    Operator curlB;

    void print_operator(const Operator &oper);

    std::vector<IndexMap> LmatX;
    std::vector<Array44> LmatX2;

    //Sources and fields on the grid
    Field3d chargeDensityOld;
    Field3d chargeDensity;
    Operator divE;
    Bounds bounds;
    

    inline int sind(int i, int j, int k) const {
        return i * ySize * zSize + j * zSize + k;
    };
    // index for 3D vector fields
    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return d + nd*(i * ySize * zSize + j * zSize + k);
    };

    inline int pos_vind(int index, int n){
        std::vector<int> dim = {xSize, ySize, zSize, 3};
        int capacity = 1;
        for(unsigned int i = n + 1; i < dim.size(); i++){
            capacity *= dim[i];
        }
        return (index / capacity) % dim[n];
    }

    void prepare();
    void computeB(const Field3d& fieldE, const Field3d& fieldEn,
                  Field3d& fieldB, double dt);
    void fdtd_explicit(Field3d& E, Field3d& B,
                       const Field3d& J, const double dt);
    void set_mirrors();
    void update_Lmat(const double3& coord, const Domain& domain, double charge,
                     double mass, double mpw, const Field3d& fieldB, const double dt);
    void update_Lmat2(const double3& coord, const Domain& domain, double charge,
                      double mass, double mpw, const Field3d& fieldB,
                      const double dt);
    void update_LmatNGP(const double3& coord, const Domain& domain, double charge,
                      double mass, double mpw, const Field3d& fieldB,
                      const double dt);



    void apply_periodic_boundaries(std::vector<IndexMap>& LmatX);
    void apply_open_boundaries_z(std::vector<IndexMap>& LmatX);
    void apply_boundaries(std::vector<IndexMap>& LmatX);

    void apply_periodic_boundaries(Field3d& field);
    void apply_open_boundaries_z(Field3d& field);
    void apply_boundaries(Field3d& field);

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

    void stencil_curlB(const Domain& domain);
    void stencil_curlE(const Domain& domain);
    void stencil_Imat(const Domain& domain);
    void stencil_Lmat(const Domain& domain);
    void stencil_divE(const Domain& domain);
    void stencil_curlE_periodic(std::vector<Trip>& trips, const Domain& domain);
    void stencil_curlE_openZ(std::vector<Trip>& trips, const Domain& domain);
    void stencil_curlB_periodic(std::vector<Trip>& trips, const Domain& domain);
    void stencil_curlB_openZ(std::vector<Trip>& trips, const Domain& domain);
    void predictE(Field3d& Ep, const Field3d& E, const Field3d& B,
                  const Field3d& J, const double dt);
    void correctE(Field3d& En, const Field3d& E, const Field3d& B,
                  const Field3d& J, const double dt);
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

    // size_t get_col_index_Lx(int i, int j, int k, int d){
        
    //     if(d == 0) return ind(i, j, k, 2, 2, 2);
    //     if (d == 1)
    //         return 8 + ind(i, j, k, 2, 3, 2);
    //     if (d == 2)
    //         return 8 + 12 + ind(i, j, k, 2, 2, 3);
    // }

    size_t get_col_index_Lx(int i, int j, int k, int d) {
        static constexpr size_t offsets[] = {
            0,       // d = 0
            8,       // d = 1
            8 + 18   // d = 2
        };

        return offsets[d] + ind2(i, j, k, 2, 2 + (d == 1), 2 + (d == 2));
    }

    size_t get_col_index_Ly(int i, int j, int k, int d) {
        static constexpr size_t offsets[] = {
            0,       // d = 0
            18,       // d = 1
            18 + 8   // d = 2
        };

        return offsets[d] +
               ind2(i, j, k, 2 + (d == 0), 2 , 2 + (d == 2));
    }

    size_t get_col_index_Lz(int i, int j, int k, int d) {
        static constexpr size_t offsets[] = {
            0,       // d = 0
            18,       // d = 1
            18 + 18   // d = 2
        };

        return offsets[d] + ind2(i, j, k, 2 + (d == 0), 2 + (d == 1), 2);
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
                if (value != 0) {
                    std::cout << ix << " " << iy << " " << iz << " " << id
                              << " " << ix1 << " " << iy1 << " " << iz1 << " "
                              << id1 << " " << value << "\n";
                }
            }
        }
    }
            void stencil_Lmat2(const Domain& domain);

           private:
            double xCellSize;
            double yCellSize;
            double zCellSize;
            int xSize;
            int ySize;
            int zSize;
        };

void print_operator(const Operator &oper);

#endif 
