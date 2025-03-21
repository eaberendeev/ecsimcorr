#ifndef COIL_H_
#define COIL_H_
#include "World.h"
struct Coil{
    double z0, R, I;
    Coil(double z0,double R,double I):z0(z0),R(R),I(I){}
};


struct CoilsArray{
    std::vector<Coil> coils;
    static const int N = 2000;
    const double hp = 2*M_PI/N;
    double R, z0, I;
    alignas(64) double cs[N];
    CoilsArray(const ParametersMap& parameters){
        auto nCoils = parameters.get_int("BCoil",0);

        for (auto k =0; k < nCoils; k++){
            z0 = parameters.get_double("BCoil",3*k + 1);
            R = parameters.get_double("BCoil", 3 * k + 2);
            I = parameters.get_double("BCoil", 3 * k + 3);

            coils.emplace_back(Coil(z0,R,I));

            std::cout << "z0 " << z0 <<", I " << I<< ", R " << R << "\n";
        }
        #pragma omp simd
        for (auto i = 0; i < N; i++){
            cs[i] = cos(i*hp);
        }
    }
    ~CoilsArray(){
    }

    double get_Bz(double z, double r);
    double get_Br(double z, double r);
    double get_integ_z(double z, double r, double R);
    double get_integ_r(double z, double r, double R);
};
void set_coils(Field3d& fieldB, const Domain& domain,
               const ParametersMap& parameters);
#endif 
