#include "timer.h"

#include <algorithm>
#include <fstream>

namespace timer {
timer globalTimer("all");

std::chrono::high_resolution_clock::time_point globalStart = std::chrono::high_resolution_clock::now();
Event events[maxEvents];
AlignedInt currEvents[maxThreads]{0};

template <typename T>
inline void putField(std::ostream& fout, const char* name, const T& val) {
    fout << '"' << name << "\": " << std::setprecision(17) << val;
}

template <typename T>
inline void putFieldString(std::ostream& fout, const char* name, const T& val) {
    fout << '"' << name << "\": \"" << val << '"';
}

void writeTimerTree(const char* filename) {
    std::ofstream fout(filename);

    fout << "[\n";

    for (int64_t thrNum = 0; thrNum < maxThreads; ++thrNum) {
        const int64_t eventsCount = currEvents[thrNum].val;

        if (thrNum != 0 && eventsCount != 0) {
            fout << ",\n";
        } else {
            fout << "\n";
        }

        for (int64_t j = 0; j < eventsCount; ++j) {
            fout << "{\n";
            const Event& event = events[thrNum * maxEventsPerThread + j];

            putFieldString(fout, "name", event.name);
            fout << ",\n";
            putFieldString(fout, "ph", "X");
            fout << ",\n";
            putField(fout, "ts", std::chrono::duration<double>(event.start - globalStart).count() * 1e6);
            fout << ",\n";
            putField(fout, "dur", std::chrono::duration<double>(event.end - event.start).count() * 1e6);
            fout << ",\n";
            putField(fout, "tid", thrNum);
            fout << ",\n";
            putField(fout, "pid", 0);

            if (j < eventsCount - 1)
                fout << "\n},\n";
            else
                fout << "}\n";
        }
    }

    fout << "]" << std::endl;
}

}   // namespace timer
