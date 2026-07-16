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

void writeFullProfile(const std::string& filename) {
    writeFullProfile(filename.c_str());
}

void writeFullProfile(const char* filename) {
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
    flatTimer timer(std::source_location::current().function_name(), size);
    return __real_malloc(size);
}

extern "C" void __wrap_free(void* p) {
    flatTimer timer(std::source_location::current().function_name());
    __real_free(p);
}

extern "C" void* __wrap_calloc(size_t n, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size * n);
    return __real_calloc(n, size);
}

extern "C" void* __wrap_realloc(void* p, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size);
    return __real_realloc(p, size);
}
extern "C" void* __wrap_reallocarray(void* p, size_t n, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size * n);
    return __real_reallocarray(p, n, size);
}

// operator new(unsigned long)
extern "C" void* __wrap__Znwm(size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size);
    return __real__Znwm(size);
}

// operator new[](unsigned long)
extern "C" void* __wrap__Znam(size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size);
    return __real__Znam(size);
}

// operator delete(void*, unsigned long)
extern "C" void __wrap__ZdlPvm(void* p, size_t size) {
    flatTimer timer(std::source_location::current().function_name(), size);
    __real__ZdlPvm(p, size);
}

#endif
#if defined(MEMORY_TIMERS_LEVEL) && MEMORY_TIMERS_LEVEL > 1

extern "C" void* __real_memcpy(void* dest, const void* src, size_t count);
extern "C" void* __real_memset(void* dest, int ch, size_t count);

extern "C" void* __wrap_memcpy(void* dest, const void* src, size_t count) {
    flatTimer timer(std::source_location::current().function_name(), count);
    return __real_memcpy(dest, src, count);
}

extern "C" void* __wrap_memset(void* dest, int ch, std::size_t count) {
    flatTimer timer(std::source_location::current().function_name(), count);
    return __real_memset(dest, ch, count);
}
#endif

}   // namespace timer
