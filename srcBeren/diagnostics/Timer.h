#ifndef TIMER_H_
#define TIMER_H_
#include <omp.h>

#include <cstdio>
#include <map>
#include <string>

#include "util.h"

struct Timer {
    FILE* fTimes;
    std::map<std::string, double> times;
    std::map<std::string, double> startT;
    bool firstWrite;
    void init(const std::string& filename) {
        create_directory(".//Performance");

        fTimes = fopen((".//Performance//" + filename).c_str(), "w");
        firstWrite = true;
    }

    void start(const std::string& s) {
        startT[s] = omp_get_wtime();
    }
    void finish(const std::string& s) {
        times[s] += omp_get_wtime() - startT[s];
    }
    void reset() {
        for (auto it = times.begin(); it != times.end(); ++it) {
            it->second = 0.;
        }
    }
    void write(int timestep, int delay);
};

#endif
