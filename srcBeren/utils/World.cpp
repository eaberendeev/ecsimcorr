#include "World.h"

Domain::Domain() {
    mCellSize = Vector3R(1.0, 1.0, 1.0);
    mOrigin = Vector3I(0, 0, 0);
    mNumCells = Vector3I(1, 1, 1);
    mSize = mNumCells + Vector3I(GHOST_NODES, GHOST_NODES, GHOST_NODES);
    bounds_.lower = {BoundType::PERIODIC, BoundType::PERIODIC,
                    BoundType::PERIODIC};
    bounds_.upper = {BoundType::PERIODIC, BoundType::PERIODIC,
                    BoundType::PERIODIC};
}

Domain::Domain(const nlohmann::json& config, const Bounds& bound) {
    set_domain(config, bound);
}

void Domain::set_domain(const nlohmann::json& config, const Bounds& bound) {
    mCellSize =
        Vector3R(get_checked<double>(config, "Dx"), get_checked<double>(config, "Dy"),
                get_checked<double>(config, "Dz"));
    mOrigin = Vector3I(0, 0, 0);
    mNumCells = Vector3I(get_checked<int>(config, "NumCellsX"),
                     get_checked<int>(config, "NumCellsY"),
                     get_checked<int>(config, "NumCellsZ"));
    mSize = mNumCells + Vector3I(GHOST_NODES, GHOST_NODES, GHOST_NODES);
    bounds_.set_bounds(config);
}

void Domain::set_domain(const nlohmann::json& config) {
    mCellSize = Vector3R(get_checked<double>(config, "Dx"),
                         get_checked<double>(config, "Dy"),
                         get_checked<double>(config, "Dz"));

    mOrigin = Vector3I(0, 0, 0);
    mNumCells = Vector3I(get_checked<int>(config, "NumCellsX"),
                         get_checked<int>(config, "NumCellsY"),
                         get_checked<int>(config, "NumCellsZ"));

    mSize = mNumCells + Vector3I(GHOST_NODES, GHOST_NODES, GHOST_NODES);
    bounds_.set_bounds(config);
}
void Domain::get_interpolation_env(const Vector3R coord, Vector3I& index, Vector3R& weight,
                           double shift = 0.0) const {
    Vector3R coordInCell =
        coord / mCellSize + Vector3R(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);

    coordInCell -= Vector3R(shift, shift, shift);

    index.x() = static_cast<int>(coordInCell.x());
    index.y() = static_cast<int>(coordInCell.y());
    index.z() = static_cast<int>(coordInCell.z());

    weight.x() = 1 - (index.x() - coordInCell.x());
    weight.y() = 1 - (index.y() - coordInCell.y());
    weight.z() = 1 - (index.z() - coordInCell.z());
}

InterpolationEnvironment Domain::get_interpolation_environment(
    const Vector3R coord, double shift = 0.0) const {
    Vector3R coordInCell =
        coord / mCellSize + Vector3R(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);

    coordInCell -= Vector3R(shift, shift, shift);

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

Vector3R Domain::interpolate_fieldB(const Field3d& field, const Vector3R& coord) {
    const int shapeSize = 2;
    Vector3R B(0,0,0);
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
