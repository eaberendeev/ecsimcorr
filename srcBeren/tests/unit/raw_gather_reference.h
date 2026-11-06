#pragma once

#include "include/config.h"
#include "utils/containers.h"

// Independent pre-migration raw kernels, retained only for regression tests.
namespace raw_gather_reference {
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

}   // namespace raw_gather_reference
