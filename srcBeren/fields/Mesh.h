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

    Field3d fieldE;
    Field3d fieldEn;
    Field3d fieldEp;
    Field3d fieldB;
    Field3d fieldJp; // predict current for EM solver
    Field3d fieldJp_full; // predict current for EM solver Jp + Lmat(E+E_n);
    Field3d fieldJe; // Esirkepov current for E correction
    Field3d fieldB0;
    Field3d fieldBInit;
    std::vector<IndexMap> LmatX;
    
    //Sources and fields on the grid
    Field3d chargeDensityOld;
    Field3d chargeDensity;
    Operator divE;

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
    void fdtd_explicit(const double dt);
    void set_mirrors();
    void update_Lmat(const double3& coord, const Domain& domain, double charge,
                     double mass, double mpw, const double dt);
    void make_periodic_border_with_add(Field3d &field );
    void make_periodic_border_with_add(Array3D<double>& field);
    void glue_Lmat_bound();

    void set_uniform_field(Field3d& field, double bx, double by, double bz);

    double3 get_fieldE_in_cell(int i, int j, int k)  const;
    double3 get_fieldB_in_cell(int i, int j, int k)  const;
    double calc_energy_field(const Field3d& field) const;
    double calc_JE(const Field3d& fieldE, const Field3d& fieldJ) const;
    double3 calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ) const;
    double get_fieldE_energy() const { return calc_energy_field(fieldE); };
    double get_fieldB_energy() const{
        return calc_energy_field(fieldB); //-calc_energy_field(fieldB0);
    };

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
    void predictE(const double dt);
    void correctE(const double dt);

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
