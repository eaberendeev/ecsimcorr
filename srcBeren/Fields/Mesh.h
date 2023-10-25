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

double3 get_fieldE_in_pos(const Field3d& fieldE, const double3& r);
double3 get_fieldB_in_pos(const Field3d& fieldB, const double3& r);
void get_fields_in_pos(const Field3d& fieldE,const Field3d& fieldB, const double3& r, double3& locE, double3 &locB);


template <typename T> 
inline void swap_add(T &a, T &b){
    a += b;
    b = a;
}



struct Mesh{
    Mesh(const World& world);

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
        return i * _size2 * _size3 + j * _size3 + k;
    };
    // index for 3D vector fields
    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return d + nd*(i * _size2 * _size3 + j * _size3 + k);
    };

    inline int pos_vind(int index, int n){
        std::vector<int> dim = {_size1, _size2, _size3, 3};
        int capacity = 1;
        for(unsigned int i = n + 1; i < dim.size(); i++){
            capacity *= dim[i];
        }
        return (index / capacity) % dim[n];
    }

    void write_field_to_file(const std::string& dataName, Field3d& field);
    void read_field_from_file(const std::string& dataName, Field3d& field);
    void write_fields_to_recovery(int timestep);
    void read_fields_from_recovery();
    void swap_add(IndexMap &a, IndexMap &b);
    void set_fields();
    void prepare();
    void computeB();
    void fdtd_explicit();
    void set_mirrors();
    void update_Lmat( const double3& coord, double charge, double mass, double mpw);
    void make_periodic_border_with_add(Field3d &field );
    void glue_Lmat_bound();

    void set_uniform_fields();

    double3 get_fieldE_in_cell(int i, int j, int k)  const;
    double3 get_fieldB_in_cell(int i, int j, int k)  const;
    double calc_energy_field(const Field3d& field) const;
    double calc_JE(const Field3d& fieldE,const Field3d& fieldJ) const;
    double get_fieldE_energy() const{
        return calc_energy_field(fieldE);
    };
    double get_fieldB_energy() const{
        return calc_energy_field(fieldB); //-calc_energy_field(fieldB0);
    };

    ~Mesh(){
    }

    double get_coord_from_nodeX(int index) const{
        return (index - CELLS_SHIFT) * Dx;
    }
    double get_coord_from_nodeY(int index) const{
        return (index - CELLS_SHIFT) * Dy;
    }
    double get_coord_from_nodeZ(int index) const{
        return (index - CELLS_SHIFT) * Dz;
    }
    int get_node_from_coordX(double coord) const{
        return int(coord / Dx + CELLS_SHIFT); 
    }
    int get_node_from_coordY(double coord) const{
        return int(coord / Dy + CELLS_SHIFT); 
    }
    int get_node_from_coordZ(double coord) const{
        return int(coord / Dz + CELLS_SHIFT); 
    }

    void stencil_curlB();
    void stencil_curlE();
    void stencil_Mmat();
    void stencil_Imat();
    void stencil_Lmat();
    void stencil_divE();
    void predictE();
    void correctE();

public:
    const World &_world;
    int _size1, _size2, _size3;
};

void print_operator(const Operator &oper);

#endif 
