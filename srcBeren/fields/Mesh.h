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
double calc_JE(const Field3d& fieldE, const Field3d& fieldJ,
               const Bounds& bounds);
double3 calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ,
                          const Bounds& bounds);
void make_periodic_border_with_add(Field3d& field, const Bounds& bounds);

// void get_fields_in_pos(const Field3d& fieldE,const Field3d& fieldB, const
// double3& r, double3& locE, double3 &locB);

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
   // void make_periodic_border_with_add(Array3D<double>& field);
    void glue_Lmat_bound();
    void set_open_bound_z(std::vector<IndexMap>& LmatX);
    void set_open_bound_z(Field3d& field);
    void set_uniform_field(
        Field3d& field, double bx, double by, double bz);

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
    // void add_init_fieldB(Field3d &B){
    //     B.data() += fieldBInit.data();
    // }
    // void remove_init_fieldB(Field3d &B){
    //     B.data() -= fieldBInit.data();
    // }

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
