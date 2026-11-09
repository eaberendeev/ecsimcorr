#ifndef COIL_H_
#define COIL_H_
#include <vector>

#include "grid/domain.h"
#include "utils/containers.h"

struct Coil {
    double z0, R, I;

    Coil(double z0, double R, double I) : z0(z0), R(R), I(I) {
    }
};

struct CoilsArray {
    std::vector<Coil> coils;
    static const int N = 2000;
    const double hp = 2 * M_PI / N;
    alignas(64) double cs[N];
    explicit CoilsArray(const std::vector<Coil>& parameters);

    double get_Bz(double z, double r);
    double get_Br(double z, double r);
    double get_integ_z(double z, double r, double R);
    double get_integ_r(double z, double r, double R);
};

void set_coils(Field3d& fieldB, const Domain& domain, const std::vector<Coil>& parameters);

#endif
