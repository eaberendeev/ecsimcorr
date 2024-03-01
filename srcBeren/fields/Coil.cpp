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

void set_coils(Field3d& fieldB, const World& world) {
    auto size_x = fieldB.size().x();
    auto size_y = fieldB.size().y();
    auto size_z = fieldB.size().z();
    double center_x = 0.5 * world.regionGlob.numCells.x() * Dx;
    double center_y = 0.5 * world.regionGlob.numCells.y() * Dy;
    double xx, yy, rr;
    double Brx, Bry, Bz;
    CoilsArray coils;

#pragma omp parallel for private(xx, yy, rr, Brx, Bry, Bz)
    for (auto k = 0; k < size_z; k++) {
        double zz = k * Dz - Dz * CELLS_SHIFT;
        for (auto i = 0; i < size_x; i++) {
            for (auto j = 0; j < size_y; j++) {


                xx = i * Dx - center_x - Dx * CELLS_SHIFT;
                yy = (j + 0.5) * Dy - center_y - Dy * CELLS_SHIFT;
                rr = sqrt(xx * xx + yy * yy);
                Brx = coils.get_Br(zz + 0.5 * Dz, rr);
                fieldB(i, j, k, 0) += Brx * xx / rr;

                yy = j * Dy - center_y - Dy * CELLS_SHIFT;
                xx = (i + 0.5) * Dx - center_x - Dx * CELLS_SHIFT;
                rr = sqrt(xx * xx + yy * yy);
                Bry = coils.get_Br(zz + 0.5 * Dz, rr);
                fieldB(i, j, k, 1) += Bry * yy / rr;
                xx = (i + 0.5) * Dx - center_x - Dx * CELLS_SHIFT;
                yy = (j + 0.5) * Dy - center_y - Dy * CELLS_SHIFT;
                rr = sqrt(xx * xx + yy * yy);
                Bz = coils.get_Bz(zz, rr);
                fieldB(i, j, k, 2) += Bz;
            }
        }
    }

    std::cout << "Coils configuration set successfull!\n";
}