#include "timer.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>

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

void writeTimerTree(const std::string& filename) {
    writeTimerTree(filename.c_str());
}

void writeTimerTree(const char* filename) {
    std::filesystem::path outFile = std::filesystem::absolute(filename);
    std::filesystem::path tmpOutFile =
        std::filesystem::absolute(filename).replace_filename("." + outFile.filename().string());
    std::ofstream fout(tmpOutFile);

    fout << "[\n";

    bool isPrintedBeforeComma = false;

    for (int64_t thrNum = 0; thrNum < maxThreads; ++thrNum) {
        const int64_t eventsCount = currEvents[thrNum].val;
        for (int64_t j = 0; j < eventsCount; ++j) {
            const Event& event = events[thrNum * maxEventsPerThread + j];

            // handle case, when some timers still unfinished before this function
            if (event.name == nullptr) {
                continue;
            }

            if (isPrintedBeforeComma) {
                fout << ",\n";
                isPrintedBeforeComma = false;
            } else {
                fout << "\n";
            }

            fout << "{\n";

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
            fout << ",\n";
            fout << "\"args\": {";
            putField(fout, "m", event.m);
            fout << "}}";

            if (!fout) {
                std::cerr << "Unable to write profile data to a temperrray file" << tmpOutFile << ", abort writing"
                          << std::endl;
                return;
            }

            isPrintedBeforeComma = true;
        }
    }

    fout << "]" << std::endl;
    fout.close();

    std::filesystem::rename(tmpOutFile, outFile);
}

void printTreeTableImpl(bool isHead, std::ostream& os, const std::vector<std::string_view>& names) {
    if (isHead) {
        for (const std::string_view& name : names) {
            const int width = std::max<int>(6, name.length());
            os << std::setw(width) << name << " ";
        }
        os << std::endl;
        return;
    }

    for (const std::string_view& name : names) {
        const double time = globalTimer.fieldTime(name);
        const int width = std::max<int>(6, name.length());
        os << std::setw(width) << std::setprecision(5) << time << " ";
    }
    os << std::endl;
}

void printSliceImpl(std::ostream& os, const std::vector<std::string_view>& names) {
    int maxNameLen = 0;
    for (const std::string_view& name : names) {
        maxNameLen = std::max<int>(maxNameLen, name.length());
    }
    for (const std::string_view& name : names) {
        const double time = globalTimer.fieldTime(name);
        const double calls = globalTimer.fieldCalls(name);
        os << std::setw(maxNameLen) << std::left << name << " " << time << "[s], " << calls << "[calls]" << std::endl;
    }
}

}   // namespace timer
