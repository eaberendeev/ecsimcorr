#ifndef COIL_H_
#define COIL_H_
#include "World.h"
struct Coil{
    double z0, R, I;
    Coil(double z0,double R,double I):z0(z0),R(R),I(I){}
};
namespace CGS
{
    static const double e   = 4.80320471257e-10;
    static const double c   = 2.99792458e10;
    static const double me  = 9.1094093837015e-28;
    static const double mc2 = 512.;
};


struct CoilsArray{
    std::vector<Coil> coils;
    const int N = 2000;
    const double hp = 2*M_PI/N;
    double R, z0, I;
    double* cs;
    CoilsArray(){
        auto nCoils = BCoil[0];
        //double n0 = 1e13;
        //double omega_p = sqrt(4.*M_PI*n0*CGS::e*CGS::e/CGS::me);
        //double B_norm = CGS::me * CGS::c * omega_p / CGS::e;
        //double B0 = 0.4 / B_norm;
        for (auto k =0; k < nCoils; k++){
            z0 = BCoil[3*k + 1];
            R = BCoil[3*k + 2];
            I = BCoil[3*k + 3];

            //double RT =  fabs(344/2) / sqrt(pow(8, 2. / 3.) - 1.);
            //R = 344/2.39;
            //double IT = 0.4 * R / (2. * M_PI);
            //R = 344/2.39; //BCoil[2+3*k];
            //I = 0.4*R/(8*M_PI)*pow(1+344*344/(R*R)/4,1.5); //BCoil[3+3*k];
            coils.emplace_back(Coil(z0,R,I));

            std::cout << "z0 " << z0 <<", I " << I<< ", R " << R << "\n";
        }
        cs = new double[N];
        for (auto i = 0; i < N; i++){
            cs[i] = cos(i*hp);
        }
    }
    ~CoilsArray(){
        delete[] cs;
    }

    double get_Bz(double z, double r);
    double get_Br(double z, double r);
    double get_integ_z(double z, double r, double R);
    double get_integ_r(double z, double r, double R);
};
void set_coils( Field3d& fieldB,const World& world);
#endif 
