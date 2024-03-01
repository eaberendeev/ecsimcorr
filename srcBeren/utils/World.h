#ifndef WORLD_H_
#define WORLD_H_
#include "Vec.h"
#include "const.h"
#include "defines.h"
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include "random_generator.h"

struct Region {
    double3 cellSize;
    int3 cellsShift;
    double origin;
    int3 numCells;
    int3 numNodes;
    int3 dampCells[2];
    int3 boundType[2];
    bool in_region(double x) const{
        if ( x < origin || x >= origin + numCells.x() * cellSize.x() )
            return false;
  
        return true;
    }
    double3 get_coord_loc(const double3 &POS) const{
        double3 POS_loc = POS;
        POS_loc.x() -= origin;
        return POS_loc;
    }
    double3 get_coord_glob(const double3 &POS) const{
        double3 POS_glob = POS;
        POS_glob.x() += origin;
        return POS_glob;    
    }
    int get_index_loc(int indx) const{
        return indx - round(origin / Dx);
    }
    Region();
    int total_size() const {
        return numNodes.x()*numNodes.y()*numNodes.z(); 
    };
};

Region split_region(const Region& regionGlob, int rank, int splitSize);

struct World{
    World(const Region& regionGlob,const Region& regionSplit): 
                                regionGlob(regionGlob),region(regionSplit){
    };
    Region regionGlob;
    Region region;
    ~World(){
    }

};

#endif 
