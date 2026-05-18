#include "timer.h"

namespace timer {
timer globalTimer("all");

std::chrono::high_resolution_clock::time_point globalStart;
Event events[maxEvents];
AlignedInt currEvents[maxThreads]{0};

}   // namespace timer
