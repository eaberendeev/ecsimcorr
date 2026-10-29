#pragma once

#include "Shape.h"
#include "config.h"
#include "containers.h"
#include "core/Vector3.h"
#include "util.h"
inline Vector3R get_fieldE_in_cell(const Field3d& fieldE, int i, int j, int k) {
    Vector3R E;
    E.x() = 0.25 * (fieldE(i, j, k, 0) + fieldE(i, j, k + 1, 0) + fieldE(i, j + 1, k, 0) + fieldE(i, j + 1, k + 1, 0));

    E.y() = 0.25 * (fieldE(i, j, k, 1) + fieldE(i, j, k + 1, 1) + fieldE(i + 1, j, k, 1) + fieldE(i + 1, j, k + 1, 1));

    E.z() = 0.125 * (fieldE(i, j, k, 2) + fieldE(i, j, k + 1, 2) + fieldE(i, j + 1, k, 2) + fieldE(i, j + 1, k + 1, 2) +
                     fieldE(i + 1, j, k, 2) + fieldE(i + 1, j, k + 1, 2) + fieldE(i + 1, j + 1, k, 2) +
                     fieldE(i + 1, j + 1, k + 1, 2));
    return E;
}

inline Vector3R get_fieldB_in_cell(const Field3d& fieldB, int i, int j, int k) {
    Vector3R B;
    B.x() = 0.25 * (fieldB(i, j, k, 0) + fieldB(i, j, k + 1, 0) + fieldB(i + 1, j, k, 0) + fieldB(i + 1, j, k + 1, 0));

    B.y() = 0.25 * (fieldB(i, j, k, 1) + fieldB(i, j, k + 1, 1) + fieldB(i, j + 1, k, 1) + fieldB(i, j + 1, k + 1, 1));

    B.z() = fieldB(i, j, k, 2);
    return B;
}

inline Vector3R interpolateE_linear(const Field3d& fieldE, const Vector3R& normalized_coord) {
    const double xx = normalized_coord.x() + GHOST_CELLS;
    const double yy = normalized_coord.y() + GHOST_CELLS;
    const double zz = normalized_coord.z() + GHOST_CELLS;

    const int indx = int(xx);
    const int indy = int(yy);
    const int indz = int(zz);

    const int indx1 = int(xx - 0.5);
    const int indy1 = int(yy - 0.5);
    const int indz1 = int(zz - 0.5);

    const double sx1 = (xx - indx);
    const double sy1 = (yy - indy);
    const double sz1 = (zz - indz);
    const double sdx1 = (xx - indx1 - 0.5);
    const double sdy1 = (yy - indy1 - 0.5);
    const double sdz1 = (zz - indz1 - 0.5);

    const double sx0 = 1. - sx1;
    const double sy0 = 1. - sy1;
    const double sz0 = 1. - sz1;
    const double sdx0 = 1. - sdx1;
    const double sdy0 = 1. - sdy1;
    const double sdz0 = 1. - sdz1;

    Vector3R E;

    E.x() =
        sdx0 * (sy0 * (sz0 * fieldE(indx1, indy, indz, 0) + sz1 * fieldE(indx1, indy, indz + 1, 0)) +
                sy1 * (sz0 * fieldE(indx1, indy + 1, indz, 0) + sz1 * fieldE(indx1, indy + 1, indz + 1, 0))) +
        sdx1 * (sy0 * (sz0 * fieldE(indx1 + 1, indy, indz, 0) + sz1 * fieldE(indx1 + 1, indy, indz + 1, 0)) +
                sy1 * (sz0 * fieldE(indx1 + 1, indy + 1, indz, 0) + sz1 * fieldE(indx1 + 1, indy + 1, indz + 1, 0)));

    E.y() =
        sx0 * (sdy0 * (sz0 * fieldE(indx, indy1, indz, 1) + sz1 * fieldE(indx, indy1, indz + 1, 1)) +
               sdy1 * (sz0 * fieldE(indx, indy1 + 1, indz, 1) + sz1 * fieldE(indx, indy1 + 1, indz + 1, 1))) +
        sx1 * (sdy0 * (sz0 * fieldE(indx + 1, indy1, indz, 1) + sz1 * fieldE(indx + 1, indy1, indz + 1, 1)) +
               sdy1 * (sz0 * fieldE(indx + 1, indy1 + 1, indz, 1) + sz1 * fieldE(indx + 1, indy1 + 1, indz + 1, 1)));

    E.z() =
        sx0 * (sy0 * (sdz0 * fieldE(indx, indy, indz1, 2) + sdz1 * fieldE(indx, indy, indz1 + 1, 2)) +
               sy1 * (sdz0 * fieldE(indx, indy + 1, indz1, 2) + sdz1 * fieldE(indx, indy + 1, indz1 + 1, 2))) +
        sx1 * (sy0 * (sdz0 * fieldE(indx + 1, indy, indz1, 2) + sdz1 * fieldE(indx + 1, indy, indz1 + 1, 2)) +
               sy1 * (sdz0 * fieldE(indx + 1, indy + 1, indz1, 2) + sdz1 * fieldE(indx + 1, indy + 1, indz1 + 1, 2)));

    return E;
}

inline Vector3R interpolateB_linear(const Field3d& fieldB, const Vector3R& normalized_coord) {
    Vector3R B;
    const double xx = normalized_coord.x() + GHOST_CELLS;
    const double yy = normalized_coord.y() + GHOST_CELLS;
    const double zz = normalized_coord.z() + GHOST_CELLS;

    const int indx = int(xx);
    const int indy = int(yy);
    const int indz = int(zz);

    const int indx1 = int(xx - 0.5);
    const int indy1 = int(yy - 0.5);
    const int indz1 = int(zz - 0.5);

    const double sx1 = (xx - indx);
    const double sy1 = (yy - indy);
    const double sz1 = (zz - indz);
    const double sdx1 = (xx - indx1 - 0.5);
    const double sdy1 = (yy - indy1 - 0.5);
    const double sdz1 = (zz - indz1 - 0.5);

    const double sx0 = 1. - sx1;
    const double sy0 = 1. - sy1;
    const double sz0 = 1. - sz1;
    const double sdx0 = 1. - sdx1;
    const double sdy0 = 1. - sdy1;
    const double sdz0 = 1. - sdz1;

    B.x() =
        sx0 * (sdy0 * (sdz0 * fieldB(indx, indy1, indz1, 0) + sdz1 * fieldB(indx, indy1, indz1 + 1, 0)) +
               sdy1 * (sdz0 * fieldB(indx, indy1 + 1, indz1, 0) + sdz1 * fieldB(indx, indy1 + 1, indz1 + 1, 0))) +
        sx1 *
            (sdy0 * (sdz0 * fieldB(indx + 1, indy1, indz1, 0) + sdz1 * fieldB(indx + 1, indy1, indz1 + 1, 0)) +
             sdy1 * (sdz0 * fieldB(indx + 1, indy1 + 1, indz1, 0) + sdz1 * fieldB(indx + 1, indy1 + 1, indz1 + 1, 0)));

    B.y() =
        sdx0 * (sy0 * (sdz0 * fieldB(indx1, indy, indz1, 1) + sdz1 * fieldB(indx1, indy, indz1 + 1, 1)) +
                sy1 * (sdz0 * fieldB(indx1, indy + 1, indz1, 1) + sdz1 * fieldB(indx1, indy + 1, indz1 + 1, 1))) +
        sdx1 *
            (sy0 * (sdz0 * fieldB(indx1 + 1, indy, indz1, 1) + sdz1 * fieldB(indx1 + 1, indy, indz1 + 1, 1)) +
             sy1 * (sdz0 * fieldB(indx1 + 1, indy + 1, indz1, 1) + sdz1 * fieldB(indx1 + 1, indy + 1, indz1 + 1, 1)));

    B.z() =
        sdx0 * (sdy0 * (sz0 * fieldB(indx1, indy1, indz, 2) + sz1 * fieldB(indx1, indy1, indz + 1, 2)) +
                sdy1 * (sz0 * fieldB(indx1, indy1 + 1, indz, 2) + sz1 * fieldB(indx1, indy1 + 1, indz + 1, 2))) +
        sdx1 * (sdy0 * (sz0 * fieldB(indx1 + 1, indy1, indz, 2) + sz1 * fieldB(indx1 + 1, indy1, indz + 1, 2)) +
                sdy1 * (sz0 * fieldB(indx1 + 1, indy1 + 1, indz, 2) + sz1 * fieldB(indx1 + 1, indy1 + 1, indz + 1, 2)));

    return B;
}

template <auto ShapeFn, int SMAX>
Vector3R interpolateE(const Field3d& fieldE, const ParticleShape<ShapeFn, SMAX>& no,
                      const ParticleShape<ShapeFn, SMAX>& sh) {
    double x_sum = 0, y_sum = 0, z_sum = 0;

#pragma omp simd collapse(3) reduction(+ : x_sum, y_sum, z_sum)
    for (int z = 0; z < SMAX; ++z) {
        for (int y = 0; y < SMAX; ++y) {
            for (int x = 0; x < SMAX; ++x) {
                const int ix = no.start_.x() + x + GHOST_CELLS;
                const int iy = no.start_.y() + y + GHOST_CELLS;
                const int iz = no.start_.z() + z + GHOST_CELLS;
                const int ix05 = sh.start_.x() + x + GHOST_CELLS;
                const int iy05 = sh.start_.y() + y + GHOST_CELLS;
                const int iz05 = sh.start_.z() + z + GHOST_CELLS;
                x_sum += fieldE(ix05, iy, iz, X) * sh(x, X) * no(y, Y) * no(z, Z);
                y_sum += fieldE(ix, iy05, iz, Y) * no(x, X) * sh(y, Y) * no(z, Z);
                z_sum += fieldE(ix, iy, iz05, Z) * no(x, X) * no(y, Y) * sh(z, Z);
            }
        }
    }
    return Vector3R(x_sum, y_sum, z_sum);
}
template <auto ShapeFn, int SMAX>
Vector3R interpolateB(const Field3d& fieldB, const ParticleShape<ShapeFn, SMAX>& no,
                      const ParticleShape<ShapeFn, SMAX>& sh) {
    double x_sum = 0, y_sum = 0, z_sum = 0;

#pragma omp simd collapse(3) reduction(+ : x_sum, y_sum, z_sum)
    for (int z = 0; z < SMAX; ++z) {
        for (int y = 0; y < SMAX; ++y) {
            for (int x = 0; x < SMAX; ++x) {
                const int ix = no.start_.x() + x + GHOST_CELLS;
                const int iy = no.start_.y() + y + GHOST_CELLS;
                const int iz = no.start_.z() + z + GHOST_CELLS;
                const int ix05 = sh.start_.x() + x + GHOST_CELLS;
                const int iy05 = sh.start_.y() + y + GHOST_CELLS;
                const int iz05 = sh.start_.z() + z + GHOST_CELLS;
                x_sum += fieldB(ix, iy05, iz05, X) * no(x, X) * sh(y, Y) * sh(z, Z);
                y_sum += fieldB(ix05, iy, iz05, Y) * sh(x, X) * no(y, Y) * sh(z, Z);
                z_sum += fieldB(ix05, iy05, iz, Z) * sh(x, X) * sh(y, Y) * no(z, Z);
            }
        }
    }
    return Vector3R(x_sum, y_sum, z_sum);
}
