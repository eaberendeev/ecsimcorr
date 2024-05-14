#include "World.h"
#include "parameters_map.h"

Domain::Domain() {
    mCellSize = double3(1.0, 1.0, 1.0);
    mOrigin = int3(0, 0, 0);
    mNumCells = int3(1, 1, 1);
    mSize = mNumCells + int3(GHOST_NODES, GHOST_NODES, GHOST_NODES);
    mBound.lowerBounds = Bounds::BoundValues(
        BoundType::PERIODIC, BoundType::PERIODIC, BoundType::PERIODIC);
    mBound.upperBounds = Bounds::BoundValues(
        BoundType::PERIODIC, BoundType::PERIODIC, BoundType::PERIODIC);
}

Domain::Domain(const ParametersMap& parameters, const Bounds& bound) {
    setDomain(parameters, bound);
}

void Domain::setDomain(const ParametersMap& parameters, const Bounds& bound) {
    mCellSize =
        double3(parameters.get_double("Dx"), parameters.get_double("Dy"),
                parameters.get_double("Dz"));
    mOrigin = int3(0, 0, 0);
    mNumCells = int3(parameters.get_int("NumCellsX_glob"),
                     parameters.get_int("NumCellsY_glob"),
                     parameters.get_int("NumCellsZ_glob"));
    mSize = mNumCells + int3(GHOST_NODES, GHOST_NODES, GHOST_NODES);
    mBound.setBounds(bound.lowerBounds, bound.upperBounds);
}

void Domain::get_interpolation_env(const double3 coord, int3& index, double3& weight,
                           double shift = 0.0) const {
    double3 coordInCell =
        coord / mCellSize + double3(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);

    coordInCell -= double3(shift, shift, shift);

    index.x() = static_cast<int>(coordInCell.x());
    index.y() = static_cast<int>(coordInCell.y());
    index.z() = static_cast<int>(coordInCell.z());

    weight.x() = 1 - (index.x() - coordInCell.x());
    weight.y() = 1 - (index.y() - coordInCell.y());
    weight.z() = 1 - (index.z() - coordInCell.z());
}

InterpolationEnvironment Domain::get_interpolation_environment(
    const double3 coord, double shift = 0.0) const {
    double3 coordInCell =
        coord / mCellSize + double3(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);

    coordInCell -= double3(shift, shift, shift);

    InterpolationEnvironment env;
    env.xIndex = static_cast<int>(coordInCell.x());
    env.yIndex = static_cast<int>(coordInCell.y());
    env.zIndex = static_cast<int>(coordInCell.z());

    env.xWeight[0] = 1 - (env.xIndex - coordInCell.x());
    env.yWeight[0] = 1 - (env.yIndex - coordInCell.y());
    env.zWeight[0] = 1 - (env.zIndex - coordInCell.z());

    env.xWeight[1] = 1 - env.xWeight[0];
    env.yWeight[1] = 1 - env.yWeight[0];
    env.zWeight[1] = 1 - env.zWeight[0];
    return env;
}

double3 Domain::interpolate_fieldB(const Field3d& field, const double3& coord) {
    const int shapeSize = 2;
    double3 B(0,0,0);
    const InterpolationEnvironment ie = get_interpolation_environment(coord);
    const InterpolationEnvironment ie05 = get_interpolation_environment(coord);

    for (int i = 0; i < shapeSize; ++i) {
        const int xIndex = ie.xIndex + i;
        const int xIndex05 = ie05.xIndex + i;
        for (int j = 0; j < shapeSize; ++j) {
            const int yIndex = ie.yIndex + j;
            const int yIndex05 = ie05.yIndex + j;
            for (int k = 0; k < shapeSize; ++k) {
                const int zIndex = ie.zIndex + k;
                const int zIndex05 = ie05.zIndex + k;
                const double wx =
                    ie.xWeight[i] * ie05.yWeight[j] * ie05.zWeight[k];
                const double wy =
                    ie05.xWeight[i] * ie.yWeight[j] * ie05.zWeight[k];
                const double wz =
                    ie05.xWeight[i] * ie05.yWeight[j] * ie.zWeight[k];
                B.x() += (wx * field(xIndex, yIndex05, zIndex05, Dim::X));
                B.y() += (wy * field(xIndex05, yIndex, zIndex05, Dim::Y));
                B.z() += (wz * field(xIndex05, yIndex05, zIndex, Dim::Z));
            }
        }
    }
    return B;
}
