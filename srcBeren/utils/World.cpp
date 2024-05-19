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

// Читаем данные из файла параметров, распознаём строки и записываем данные в Params
Region::Region(){

    cellSize = double3(Dx,Dy,Dz);
    cellsShift = int3(CELLS_SHIFT,CELLS_SHIFT,CELLS_SHIFT);
    origin = 0.0;
   
    numCells = int3(NumCellsX_glob, NumCellsY_glob, NumCellsZ_glob);    
    numNodes = int3(numCells.x() + ADD_NODES,
                     numCells.y() + ADD_NODES,
                     numCells.z() + ADD_NODES);
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
