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

// Читаем данные из файла параметров, распознаём строки и записываем данные в Params
Region::Region(const ParametersMap& parameters) {
    cellSize = double3(parameters.get_double("Dx"), parameters.get_double("Dy"),
                       parameters.get_double("Dz"));
    cellsShift = int3(CELLS_SHIFT,CELLS_SHIFT,CELLS_SHIFT);
    origin = 0.0;

    numCells = int3(parameters.get_int("NumCellsX_glob"),
                    parameters.get_int("NumCellsY_glob"),
                    parameters.get_int("NumCellsZ_glob"));
    numNodes = int3(numCells.x() + GHOST_NODES, numCells.y() + GHOST_NODES,
                    numCells.z() + GHOST_NODES);
    dampCells[0] = int3(DampCellsX_glob[0],DampCellsY_glob[0],DampCellsZ_glob[0]);
    dampCells[1] = int3(DampCellsX_glob[1],DampCellsY_glob[1],DampCellsZ_glob[1]);

    boundType[0] = int3(BoundTypeX_glob[0],BoundTypeY_glob[0],BoundTypeZ_glob[0]);
    boundType[1] = int3(BoundTypeX_glob[1],BoundTypeY_glob[1],BoundTypeZ_glob[1]);
}

Region split_region(const Region& regionGlob, int rank, int splitSize){
    Region region = regionGlob;
    
    std::cout << "Current number of cells = "<< region.numCells.x() << ". Start coord =" << region.origin << std::endl;
    return region;
}
