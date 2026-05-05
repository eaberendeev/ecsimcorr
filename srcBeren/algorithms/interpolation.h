#pragma once

#include "containers.h"
#include "vector3.h"
#include "Shape.h"
#include "config.h"
#include "util.h"
inline Vector3R get_fieldE_in_cell(const Field3d& fieldE, int i, int j,
                                 int k) {
    Vector3R E;
    E.x() = 0.25 * (fieldE(i, j, k, 0) + fieldE(i, j, k + 1, 0) +
                    fieldE(i, j + 1, k, 0) + fieldE(i, j + 1, k + 1, 0));

    E.y() = 0.25 * (fieldE(i, j, k, 1) + fieldE(i, j, k + 1, 1) +
                    fieldE(i + 1, j, k, 1) + fieldE(i + 1, j, k + 1, 1));

    E.z() =
        0.125 * (fieldE(i, j, k, 2) + fieldE(i, j, k + 1, 2) +
                 fieldE(i, j + 1, k, 2) + fieldE(i, j + 1, k + 1, 2) +
                 fieldE(i + 1, j, k, 2) + fieldE(i + 1, j, k + 1, 2) +
                 fieldE(i + 1, j + 1, k, 2) + fieldE(i + 1, j + 1, k + 1, 2));
    return E;
}

inline Vector3R get_fieldB_in_cell(const Field3d& fieldB, int i, int j,
                                 int k) {
    Vector3R B;
    B.x() = 0.25 * (fieldB(i, j, k, 0) + fieldB(i, j, k + 1, 0) +
                    fieldB(i + 1, j, k, 0) + fieldB(i + 1, j, k + 1, 0));

    B.y() = 0.25 * (fieldB(i, j, k, 1) + fieldB(i, j, k + 1, 1) +
                    fieldB(i, j + 1, k, 1) + fieldB(i, j + 1, k + 1, 1));

    B.z() = fieldB(i, j, k, 2);
    return B;
}

inline Vector3R interpolateE_ngp(const Field3d& fieldE,
                                const Vector3R& normalized_coord) {
    const auto x05 = normalized_coord.x() - 0.5;
    const auto y05 = normalized_coord.y() - 0.5;
    const auto z05 = normalized_coord.z() - 0.5;

    const auto ix = ngp(normalized_coord.x() + GHOST_CELLS);
    const auto iy = ngp(normalized_coord.y() + GHOST_CELLS);
    const auto iz = ngp(normalized_coord.z() + GHOST_CELLS);
    const auto ix05 = ngp(x05 + GHOST_CELLS);
    const auto iy05 = ngp(y05 + GHOST_CELLS);
    const auto iz05 = ngp(z05 + GHOST_CELLS);

    Vector3R E;
    E.x() = fieldE(ix05, iy, iz, X);
    E.y() = fieldE(ix, iy05, iz, Y);
    E.z() = fieldE(ix, iy, iz05, Z);
    return E;
}

inline Vector3R interpolateB_ngp(const Field3d& fieldB,
                                const Vector3R& normalized_coord) {
    const auto x05 = normalized_coord.x() - 0.5;
    const auto y05 = normalized_coord.y() - 0.5;
    const auto z05 = normalized_coord.z() - 0.5;

    const auto ix = ngp(normalized_coord.x() + GHOST_CELLS);
    const auto iy = ngp(normalized_coord.y() + GHOST_CELLS);
    const auto iz = ngp(normalized_coord.z() + GHOST_CELLS);
    const auto ix05 = ngp(x05 + GHOST_CELLS);
    const auto iy05 = ngp(y05 + GHOST_CELLS);
    const auto iz05 = ngp(z05 + GHOST_CELLS);

    Vector3R B;
    B.x() = fieldB(ix, iy05, iz05, X);
    B.y() = fieldB(ix05, iy, iz05, Y);
    B.z() = fieldB(ix05, iy05, iz, Z);
    return B;
}

inline Vector3R interpolateE_linear(const Field3d& fieldE,
                                   const Vector3R& normalized_coord) {
    Vector3R E;

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

    E.x() = sdx0 * (sy0 * (sz0 * fieldE(indx1, indy, indz, 0) +
                           sz1 * fieldE(indx1, indy, indz + 1, 0)) +
                    sy1 * (sz0 * fieldE(indx1, indy + 1, indz, 0) +
                           sz1 * fieldE(indx1, indy + 1, indz + 1, 0))) +
            sdx1 * (sy0 * (sz0 * fieldE(indx1 + 1, indy, indz, 0) +
                           sz1 * fieldE(indx1 + 1, indy, indz + 1, 0)) +
                    sy1 * (sz0 * fieldE(indx1 + 1, indy + 1, indz, 0) +
                           sz1 * fieldE(indx1 + 1, indy + 1, indz + 1, 0)));

    E.y() = sx0 * (sdy0 * (sz0 * fieldE(indx, indy1, indz, 1) +
                           sz1 * fieldE(indx, indy1, indz + 1, 1)) +
                   sdy1 * (sz0 * fieldE(indx, indy1 + 1, indz, 1) +
                           sz1 * fieldE(indx, indy1 + 1, indz + 1, 1))) +
            sx1 * (sdy0 * (sz0 * fieldE(indx + 1, indy1, indz, 1) +
                           sz1 * fieldE(indx + 1, indy1, indz + 1, 1)) +
                   sdy1 * (sz0 * fieldE(indx + 1, indy1 + 1, indz, 1) +
                           sz1 * fieldE(indx + 1, indy1 + 1, indz + 1, 1)));

    E.z() = sx0 * (sy0 * (sdz0 * fieldE(indx, indy, indz1, 2) +
                          sdz1 * fieldE(indx, indy, indz1 + 1, 2)) +
                   sy1 * (sdz0 * fieldE(indx, indy + 1, indz1, 2) +
                          sdz1 * fieldE(indx, indy + 1, indz1 + 1, 2))) +
            sx1 * (sy0 * (sdz0 * fieldE(indx + 1, indy, indz1, 2) +
                          sdz1 * fieldE(indx + 1, indy, indz1 + 1, 2)) +
                   sy1 * (sdz0 * fieldE(indx + 1, indy + 1, indz1, 2) +
                          sdz1 * fieldE(indx + 1, indy + 1, indz1 + 1, 2)));

    return E;
}

inline Vector3R interpolateB_linear(const Field3d& fieldB,
                                   const Vector3R& normalized_coord) {
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

    B.x() = sx0 * (sdy0 * (sdz0 * fieldB(indx, indy1, indz1, 0) +
                           sdz1 * fieldB(indx, indy1, indz1 + 1, 0)) +
                   sdy1 * (sdz0 * fieldB(indx, indy1 + 1, indz1, 0) +
                           sdz1 * fieldB(indx, indy1 + 1, indz1 + 1, 0))) +
            sx1 * (sdy0 * (sdz0 * fieldB(indx + 1, indy1, indz1, 0) +
                           sdz1 * fieldB(indx + 1, indy1, indz1 + 1, 0)) +
                   sdy1 * (sdz0 * fieldB(indx + 1, indy1 + 1, indz1, 0) +
                           sdz1 * fieldB(indx + 1, indy1 + 1, indz1 + 1, 0)));

    B.y() = sdx0 * (sy0 * (sdz0 * fieldB(indx1, indy, indz1, 1) +
                           sdz1 * fieldB(indx1, indy, indz1 + 1, 1)) +
                    sy1 * (sdz0 * fieldB(indx1, indy + 1, indz1, 1) +
                           sdz1 * fieldB(indx1, indy + 1, indz1 + 1, 1))) +
            sdx1 * (sy0 * (sdz0 * fieldB(indx1 + 1, indy, indz1, 1) +
                           sdz1 * fieldB(indx1 + 1, indy, indz1 + 1, 1)) +
                    sy1 * (sdz0 * fieldB(indx1 + 1, indy + 1, indz1, 1) +
                           sdz1 * fieldB(indx1 + 1, indy + 1, indz1 + 1, 1)));

    B.z() = sdx0 * (sdy0 * (sz0 * fieldB(indx1, indy1, indz, 2) +
                            sz1 * fieldB(indx1, indy1, indz + 1, 2)) +
                    sdy1 * (sz0 * fieldB(indx1, indy1 + 1, indz, 2) +
                            sz1 * fieldB(indx1, indy1 + 1, indz + 1, 2))) +
            sdx1 * (sdy0 * (sz0 * fieldB(indx1 + 1, indy1, indz, 2) +
                            sz1 * fieldB(indx1 + 1, indy1, indz + 1, 2)) +
                    sdy1 * (sz0 * fieldB(indx1 + 1, indy1 + 1, indz, 2) +
                            sz1 * fieldB(indx1 + 1, indy1 + 1, indz + 1, 2)));

    return B;
}

inline Vector3R interpolateE(const Field3d& fieldE,
                            const Vector3R& normalized_coord, ShapeType type) {
    switch (type) {
        case ShapeType::NGP:
            return interpolateE_ngp(fieldE, normalized_coord);
        case ShapeType::Linear:
            return interpolateE_linear(fieldE, normalized_coord);
        case ShapeType::Quadratic:
            std::cout << "interpolateE for quadratic is not supported\n"
                      << std::endl;
    }
    return Vector3R(0, 0, 0);
}

template <auto ShapeFn, int SMAX>
Vector3R interpolateE(const Field3d& fieldE, ParticleShape<ShapeFn, SMAX>& no,
                      ParticleShape<ShapeFn, SMAX>& sh) {
    Vector3R E = Vector3R(0, 0, 0);

#pragma omp simd collapse(3)
    for (int z = 0; z < SMAX; ++z) {
        for (int y = 0; y < SMAX; ++y) {
            for (int x = 0; x < SMAX; ++x) {
                const int ix = no.start.x() + x + GHOST_CELLS;
                const int iy = no.start.y() + y + GHOST_CELLS;
                const int iz = no.start.z() + z + GHOST_CELLS;
                const int ix05 = sh.start.x() + x + GHOST_CELLS;
                const int iy05 = sh.start.y() + y + GHOST_CELLS;
                const int iz05 = sh.start.z() + z + GHOST_CELLS;
                E.x() +=
                    fieldE(ix05, iy, iz, X) * sh(x, X) * no(y, Y) * no(z, Z);
                E.y() +=
                    fieldE(ix, iy05, iz, Y) * no(x, X) * sh(y, Y) * no(z, Z);
                E.z() +=
                    fieldE(ix, iy, iz05, Z) * no(x, X) * no(y, Y) * sh(z, Z);
            }
        }
    }
    return E;
}
template <auto ShapeFn, int SMAX>
Vector3R interpolateB(const Field3d& fieldB, ParticleShape<ShapeFn, SMAX>& no,
                      ParticleShape<ShapeFn, SMAX>& sh) {
    Vector3R B = Vector3R(0, 0, 0);

#pragma omp simd collapse(3)
    for (int z = 0; z < SMAX; ++z) {
        for (int y = 0; y < SMAX; ++y) {
            for (int x = 0; x < SMAX; ++x) {
                const int ix = no.start_.x() + x + GHOST_CELLS;
                const int iy = no.start_.y() + y + GHOST_CELLS;
                const int iz = no.start_.z() + z + GHOST_CELLS;
                const int ix05 = sh.start_.x() + x + GHOST_CELLS;
                const int iy05 = sh.start_.y() + y + GHOST_CELLS;
                const int iz05 = sh.start_.z() + z + GHOST_CELLS;
                B.x() +=
                    fieldB(ix, iy05, iz05, X) * no(x, X) * sh(y, Y) * sh(z, Z);
                B.y() +=
                    fieldB(ix05, iy, iz05, Y) * sh(x, X) * no(y, Y) * sh(z, Z);
                B.z() +=
                    fieldB(ix05, iy05, iz, Z) * sh(x, X) * sh(y, Y) * no(z, Z);
            }
        }
    }
    return B;
}
