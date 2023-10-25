#include "World.h"

#define RANDOMSTD

#ifdef RANDOMSTD
  #include <random>

  std::mt19937 gen;
  std::uniform_real_distribution<> urd(0, 1); 

  double Uniform01(){
    return urd(gen);
  }

  void SetRandSeed(int val){
    gen.seed(val);
  }
#else
  double Uniform01(){
    return (double)(rand())/RAND_MAX ;
  }
  void SetRandSeed(int val){
      srand(val);
  }
#endif

double Gauss(double sigma){
  double r1 = Uniform01();
  double r2 = Uniform01();
  
  return sigma*sqrt(-2.0*log(r1))*sin(2.0*PI*r2);
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
