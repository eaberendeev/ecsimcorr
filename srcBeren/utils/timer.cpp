#include "timer.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>

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

void writeFullProfile(const std::string& filename) {
    writeFullProfile(filename.c_str());
}

struct JsonAuxPrinter {
    JsonAuxPrinter(const std::filesystem::path& path) : fout(path) {
        fout << "[\n";
        isBlockEmpty.push_back(true);
    }
    ~JsonAuxPrinter() {
        assert(isBlockEmpty.size() == 1);
        fout << "]";
    }

    inline void startBlock() {
        if (isBlockEmpty.back() == false) {
            fout << ",\n";
        }
        isBlockEmpty.back() = false;
        fout << "{\n";
        isBlockEmpty.push_back(true);
    }

    inline void startBlockNamed(const char* name) {
        if (isBlockEmpty.back() == false) {
            fout << ",\n";
        }
        isBlockEmpty.back() = false;
        fout << '"' << name << "\": " << "{\n";
        isBlockEmpty.push_back(true);
    }

    inline void finishBlock() {
        fout << "\n}";
        isBlockEmpty.pop_back();
    }

    template <typename T>
    inline void putField(const char* name, const T& val) {
        if (!isBlockEmpty.back()) {
            fout << ",\n";
        }
        fout << '"' << name << "\": " << std::setprecision(17) << val;
        isBlockEmpty.back() = false;
    }

    inline void putTimeNanoSeconds(const char* name, int64_t val) {
        if (!isBlockEmpty.back()) {
            fout << ",\n";
        }
        const char oldFill = fout.fill();

        fout << '"' << name << "\": ";
        if (val < 0) {
            fout << '-';
        }
        val = std::abs(val);
        fout << val / 1000 << "." << std::setw(3) << std::setfill('0') << val % 1000 << std::setfill(oldFill);
        isBlockEmpty.back() = false;
    }

    template <typename T>
    inline void putFieldQuoted(const char* name, const T& val) {
        if (!isBlockEmpty.back()) {
            fout << ",\n";
        }
        fout << '"' << name << "\": \"" << val << '"';
        isBlockEmpty.back() = false;
    }

    std::ofstream fout;
    std::vector<bool> isBlockEmpty;
};

void writeFullProfile(const char* filename) {
    std::filesystem::path outFile = std::filesystem::absolute(filename);
    std::filesystem::path tmpOutFile =
        std::filesystem::absolute(filename).replace_filename("." + outFile.filename().string());

    JsonAuxPrinter jsonPrinter(tmpOutFile);

    std::vector<double> prevBandwidths({0.0});
    std::vector<double> prevEndTimes({std::numeric_limits<double>::max()});

    for (int64_t thrNum = 0; thrNum < maxThreads; ++thrNum) {
        const std::string nameBandwidthOmp = "bandwidth for thr" + std::to_string(thrNum);
        const std::string nameBandwidthMaster = "bandwidth for master thr";
        const int64_t eventsCount = currEvents[thrNum].val;

        for (int64_t j = 0; j < eventsCount; ++j) {
            const Event& event = events[thrNum * maxEventsPerThread + j];

            // handle case, when some timers still unfinished before this function
            if (event.name == nullptr) {
                continue;
            }

            using NanoSecDuration = std::chrono::duration<double, std::ratio<1L, 1'000'000'000L>>;
            const int64_t start = NanoSecDuration(event.start - globalStart).count();
            const int64_t duration = NanoSecDuration(event.end - event.start).count();
            const int64_t end = NanoSecDuration(event.end - globalStart).count();

            const double gb = event.unit == MeasureUnit::byte ? event.m / 1024.0 / 1024.0 / 1024.0 : 0.0;
            const double bandwidth = gb / (1e9 * duration);

            jsonPrinter.startBlock();
            jsonPrinter.putFieldQuoted("name", event.name);
            jsonPrinter.putFieldQuoted("ph", 'X');
            jsonPrinter.putTimeNanoSeconds("ts", start);
            jsonPrinter.putTimeNanoSeconds("dur", duration);
            jsonPrinter.putField("tid", thrNum);
            jsonPrinter.putField("pid", 0);
            jsonPrinter.startBlockNamed("args");
            if (event.unit == MeasureUnit::byte) {
                jsonPrinter.putField("size Gb", gb);
                if (event.m != -1) {
                    jsonPrinter.putField("bandwidth Gb/s ", bandwidth);
                }
            } else {
                jsonPrinter.putField("m", event.m);
                if (event.m != -1) {
                    jsonPrinter.putField("m", event.m);
                    jsonPrinter.putField("perf", static_cast<double>(event.m) / duration);
                }
            }
            jsonPrinter.finishBlock();
            jsonPrinter.finishBlock();

            if (thrNum == 0 && event.isOmp) {
                continue;
            }

            assert(thrNum == 0 || event.isOmp);
            const std::string& nameBandwidth = event.isOmp ? nameBandwidthOmp : nameBandwidthMaster;

            jsonPrinter.startBlock();
            jsonPrinter.putFieldQuoted("name", nameBandwidth);
            jsonPrinter.putFieldQuoted("ph", 'C');
            jsonPrinter.putTimeNanoSeconds("ts", start);
            jsonPrinter.putField("tid", thrNum);
            jsonPrinter.putField("pid", 0);
            jsonPrinter.startBlockNamed("args");
            jsonPrinter.putField("value", bandwidth);
            jsonPrinter.finishBlock();
            jsonPrinter.finishBlock();

            while (end > prevEndTimes.back()) {
                prevBandwidths.pop_back();
                prevEndTimes.pop_back();
            }

            jsonPrinter.startBlock();
            jsonPrinter.putFieldQuoted("name", nameBandwidth);
            jsonPrinter.putFieldQuoted("ph", 'C');
            jsonPrinter.putTimeNanoSeconds("ts", end);
            jsonPrinter.putField("tid", thrNum);
            jsonPrinter.putField("pid", 0);
            jsonPrinter.startBlockNamed("args");
            jsonPrinter.putField("value", prevBandwidths.back());
            jsonPrinter.finishBlock();
            jsonPrinter.finishBlock();

            prevBandwidths.push_back(bandwidth);
            prevEndTimes.push_back(end);

            if (!jsonPrinter.fout) {
                std::cerr << "Unable to write profile data to a temperrray file" << tmpOutFile << ", abort writing"
                          << std::endl;
                return;
            }
        }
    }

    std::filesystem::rename(tmpOutFile, outFile);
}

void clearFullProfile() {
    for (int64_t thrNum = 0; thrNum < maxThreads; ++thrNum) {
        const int64_t eventsCount = currEvents[thrNum].val;
        for (int64_t j = 0; j < eventsCount; ++j) {
            events[thrNum * maxEventsPerThread + j] = Event{};
        }
        currEvents[thrNum].val = 0;
    }
}

void printTableFromTreeImpl(bool isHead, std::ostream& os, const std::vector<std::string_view>& names) {
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

#if defined(MEMORY_TIMERS_LEVEL) && MEMORY_TIMERS_LEVEL > 0

extern "C" void* __real_malloc(size_t size);
extern "C" void __real_free(void* p);
extern "C" void* __real_calloc(size_t n, size_t size);
extern "C" void* __real_realloc(void* p, size_t size);
extern "C" void* __real_reallocarray(void* p, size_t n, size_t size);
// operator new(unsigned long)
extern "C" void* __real__Znwm(size_t size);
// operator new[](unsigned long)
extern "C" void* __real__Znam(size_t size);
// operator delete(void*, unsigned long)
extern "C" void __real__ZdlPvm(void* p, size_t size);

extern "C" void* __wrap_malloc(size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size, MeasureUnit::byte);
    return __real_malloc(size);
}

extern "C" void __wrap_free(void* p) {
    flatTimer timer(std::source_location::current().function_name());
    __real_free(p);
}

extern "C" void* __wrap_calloc(size_t n, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size * n, MeasureUnit::byte);
    return __real_calloc(n, size);
}

extern "C" void* __wrap_realloc(void* p, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size, MeasureUnit::byte);
    return __real_realloc(p, size);
}
extern "C" void* __wrap_reallocarray(void* p, size_t n, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size * n, MeasureUnit::byte);
    return __real_reallocarray(p, n, size);
}

// operator new(unsigned long)
extern "C" void* __wrap__Znwm(size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size, MeasureUnit::byte);
    return __real__Znwm(size);
}

// operator new[](unsigned long)
extern "C" void* __wrap__Znam(size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size, MeasureUnit::byte);
    return __real__Znam(size);
}

// operator delete(void*, unsigned long)
extern "C" void __wrap__ZdlPvm(void* p, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size, MeasureUnit::byte);
    __real__ZdlPvm(p, size);
}

#endif
#if defined(MEMORY_TIMERS_LEVEL) && MEMORY_TIMERS_LEVEL > 1

extern "C" void* __real_memcpy(void* dest, const void* src, size_t count);
extern "C" void* __real_memset(void* dest, int ch, size_t count);

extern "C" void* __wrap_memcpy(void* dest, const void* src, size_t count) {
    flatTimer timer(std::source_location::current().function_name(), count, MeasureUnit::byte);
    return __real_memcpy(dest, src, count);
}

extern "C" void* __wrap_memset(void* dest, int ch, std::size_t count) {
    flatTimer timer(std::source_location::current().function_name(), count, MeasureUnit::byte);
    return __real_memset(dest, ch, count);
}
#endif

}   // namespace timer
