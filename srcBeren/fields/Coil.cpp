#include "Coil.h"

#include "runtime/profiler.h"

#define EPS 1.e-10

CoilsArray::CoilsArray(const std::vector<Coil>& parameters) : coils(parameters) {
#pragma omp simd
    for (int i = 0; i < N; ++i) {
        cs[i] = cos(i * hp);
    }
}

double CoilsArray::get_integ_r(double z, double r, double R) {
    RECORD_TIMER;

    double sum = 0;
    double znam;
#pragma omp simd reduction(+ : sum)
    for (auto i = 0; i < N; i++) {
        znam = R * R + z * z + r * r - 2 * R * r * cs[i];
        znam = (fabs(znam) < EPS) ? EPS : znam;
        const double denom = znam * sqrt(znam);

        sum += hp * (cs[i] / denom);
    }
    return sum;
}
double CoilsArray::get_integ_z(double z, double r, double R) {
    double sum = 0;
    double znam;

#pragma omp simd reduction(+ : sum)
    for (auto i = 0; i < N; i++) {
        znam = R * R + z * z + r * r - 2 * R * r * cs[i];
        znam = (fabs(znam) < EPS) ? EPS : znam;
        const double denom = znam * sqrt(znam);

        sum += hp * ((R - r * cs[i]) / denom);
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

void set_coils(Field3d& fieldB, const Domain& domain, const std::vector<Coil>& parameters) {
    CoilsArray coils(parameters);
    auto size_x = fieldB.sizes().x();
    auto size_y = fieldB.sizes().y();
    auto size_z = fieldB.sizes().z();
    const Vector3R center = 0.5 * (domain.geom.box_min + domain.geom.box_max);
    const double center_x = center.x(), center_y = center.y();
    double xx, yy, rr;
    double Brx, Bry, Bz;

#pragma omp parallel for collapse(3) private(xx, yy, rr, Brx, Bry, Bz)
    for (auto k = 0; k < size_z; k++) {
        for (auto i = 0; i < size_x; i++) {
            for (auto j = 0; j < size_y; j++) {
                const Vector3R bx_node = domain.get_node_position(i, j, k, FieldType::MAGNETIC, 0);
                xx = bx_node.x() - center_x;
                yy = bx_node.y() - center_y;
                rr = std::hypot(xx, yy);
                Brx = coils.get_Br(bx_node.z(), rr);
                // Bx  = cos(phi)*Br
                fieldB(i, j, k, 0) += rr == 0 ? 0 : Brx * xx / rr;

                const Vector3R by_node = domain.get_node_position(i, j, k, FieldType::MAGNETIC, 1);
                yy = by_node.y() - center_y;
                xx = by_node.x() - center_x;
                rr = std::hypot(xx, yy);
                Bry = coils.get_Br(by_node.z(), rr);
                // By  = sin(phi)*Br
                fieldB(i, j, k, 1) += rr == 0 ? 0 : Bry * yy / rr;
                const Vector3R bz_node = domain.get_node_position(i, j, k, FieldType::MAGNETIC, 2);
                xx = bz_node.x() - center_x;
                yy = bz_node.y() - center_y;
                rr = std::hypot(xx, yy);
                Bz = coils.get_Bz(bz_node.z(), rr);
                fieldB(i, j, k, 2) += Bz;
            }
        }
    }

    std::cout << "Coils configuration set successfull!\n";
}
