#include "Coil.h"

double CoilsArray::get_integ_r(double z, double r, double R) {
    double sum = 0;
    double znam;
    for (auto i = 0; i < N; i++) {
        znam = R * R + z * z + r * r - 2 * R * r * cs[i];
        sum += hp * (cs[i] / (znam * sqrt(znam)));
    }
    return sum;
}
double CoilsArray::get_integ_z(double z, double r, double R) {
    double sum = 0;
    double znam;

    for (auto i = 0; i < N; i++) {
        znam = R * R + z * z + r * r - 2 * R * r * cs[i];
        sum += hp * ((R - r * cs[i]) / (znam * sqrt(znam)));
    }
    return sum;
}

double CoilsArray::get_Bz(double z, double r) {
    double Bz = 0;
    for (const auto& coil : coils) {
        double zc = z - coil.z0;
        Bz += coil.I * coil.R * get_integ_z(zc, r, coil.R);
    }
    return Bz;
}
double CoilsArray::get_Br(double z, double r) {
    double Br = 0;
    for (const auto& coil : coils) {
        double zc = z - coil.z0;
        Br += coil.I * coil.R * zc * get_integ_r(zc, r, coil.R);
    }
    return Br;
}

void set_coils(Field3d& fieldB, const Domain &domain, const ParametersMap &parameters) {
    auto size_x = fieldB.size().x();
    auto size_y = fieldB.size().y();
    auto size_z = fieldB.size().z();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();
    double center_x = 0.5 * size_x * dx;
    double center_y = 0.5 * size_y * dy;
    double xx, yy, rr;
    double Brx, Bry, Bz;
    CoilsArray coils(parameters);

#pragma omp parallel for private(xx, yy, rr, Brx, Bry, Bz)
    for (auto k = 0; k < size_z; k++) {
        double zz = k * dz - dz * GHOST_CELLS;
        for (auto i = 0; i < size_x; i++) {
            for (auto j = 0; j < size_y; j++) {

// TO DO: Get cordinate of field in nodes
                xx = i * dx - center_x - dx * GHOST_CELLS;
                yy = (j + 0.5) * dy - center_y - dy * GHOST_CELLS;
                rr = sqrt(xx * xx + yy * yy);
                Brx = coils.get_Br(zz + 0.5 * dz, rr);
                fieldB(i, j, k, 0) += Brx * xx / rr;

                yy = j * dy - center_y - dy * GHOST_CELLS;
                xx = (i + 0.5) * dx - center_x - dx * GHOST_CELLS;
                rr = sqrt(xx * xx + yy * yy);
                Bry = coils.get_Br(zz + 0.5 * dz, rr);
                fieldB(i, j, k, 1) += Bry * yy / rr;
                xx = (i + 0.5) * dx - center_x - dx * GHOST_CELLS;
                yy = (j + 0.5) * dy - center_y - dy * GHOST_CELLS;
                rr = sqrt(xx * xx + yy * yy);
                Bz = coils.get_Bz(zz, rr);
                fieldB(i, j, k, 2) += Bz;
            }
        }
    }

    std::cout << "Coils configuration set successfull!\n";
}
