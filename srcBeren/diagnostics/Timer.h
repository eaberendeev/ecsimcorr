#ifndef TIMER_H_
#define TIMER_H_
#include "World.h"
#include <omp.h>

struct Timer{
    FILE *fTimes;
    std::map<std::string,double> times;
    std::map<std::string,double> startT;
    bool firstWrite;
    Timer(const std::string& filename){
        fTimes = fopen( (".//Performance//"+filename).c_str(), "w");
        _reset = true;
        firstWrite = true;
    }

    void start(const std::string &s){
        startT[s] = omp_get_wtime();
    }
    void finish(const std::string &s){

        times[s] += omp_get_wtime() - startT[s];
    }
    void reset(){
        for (auto it = times.begin(); it != times.end(); ++it){
            it->second = 0.;
        }
        _reset = true;
    }
    void write(int timestep, int delay); 
protected:
    bool _reset;
};

#endif 	
