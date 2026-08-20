#pragma once

#include <omp.h>

#include <chrono>
#include <cstdint>
#include <iostream>

namespace timer {

enum class MeasureUnit : uint8_t {
    not_set,
    dim,
    byte,
};

struct NoStart {};

}   // namespace timer

#ifdef USE_TIMERS

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <source_location>
#include <string_view>
#include <vector>

namespace timer {

class timer {
   public:
    timer(const char* name) : name(name), start(now()) {
        if (!omp_in_parallel())
            init();
    }

    void finish() {
        if (!isInitialized)
            return;

        if (isFinished)
            return;

        isFinished = true;
        end = now();
        duration = std::chrono::duration<double>(this->end - this->start).count();

        if (prevTimer != nullptr) {
            activeTimer() = prevTimer;
            std::vector<timer>::iterator begin = prevTimer->lowerTimers.begin();
            std::vector<timer>::iterator end = prevTimer->lowerTimers.end();
            std::vector<timer>::iterator pos = std::find(begin, end, *this);
            calls = 1;
            if (pos == end)
                prevTimer->lowerTimers.push_back(*this);
            else
                pos->mergeTimers(*this);
        }
    }

    ~timer() {
        finish();
    }

    bool operator==(const timer& over) const {
        return name == over.name || std::string_view(name) == over.name;
    }

    void printTimers(const int64_t nestingDepth, std::ostream& os) const {
        for (int64_t i = 0; i < nestingDepth; ++i) {
            os << "|  ";
        }
        const std::string_view cuttedName = cutName();
        os << "> " << cuttedName << ": " << calls << "[calls], " << duration << "[s]";
        if (calls != 1) {
            os << ", " << duration / calls << "[s/call]";
        }
        os << std::endl;
        double accumulator = 0.0;
        for (const auto& it : lowerTimers) {
            it.printTimers(nestingDepth + 1, os);
            accumulator += it.duration;
        }
        if (std::ssize(lowerTimers)) {
            for (int64_t i = 0; i < nestingDepth; ++i) os << "|  ";
            os << "\\  > " << "self " << ": " << duration - accumulator << "[s]" << std::endl;
        }
    }

    double fieldTime(const std::string_view reqName) const {
        double sum = 0.0;
        // std::cout<<"cutName : "<<cutName()<<" req name: "<<reqName<<std::endl;
        if (cutName() == reqName)
            sum += duration;
        for (const timer& lower : lowerTimers) {
            sum += lower.fieldTime(reqName);
        }

        return sum;
    }

    int fieldCalls(const std::string_view reqName) const {
        int sum = 0.0;
        // std::cout<<"cutName : "<<cutName()<<" req name: "<<reqName<<std::endl;
        if (cutName() == reqName)
            sum += calls;
        for (const timer& lower : lowerTimers) {
            sum += lower.fieldCalls(reqName);
        }

        return sum;
    }

    void clear() {
        isFinished = false;
        duration = 0.0;
        calls = 1;
        lowerTimers.clear();
        start = now();
    }

   private:
    std::string_view cutName() const {
        const std::string_view nameView = name;
        const size_t endPos = nameView.find('(');
        const std::string_view cuttedName = endPos == std::string_view::npos ? nameView : nameView.substr(0, endPos);
        return cuttedName;
    }

    void mergeTimers(const timer& over) {
        assert(*this == over);
        duration += over.duration;
        calls += over.calls;

        for (const auto& it : over.lowerTimers) {
            std::vector<timer>::iterator begin = lowerTimers.begin();
            std::vector<timer>::iterator end = lowerTimers.end();
            std::vector<timer>::iterator pos = std::find(begin, end, it);
            if (pos == end)
                lowerTimers.push_back(it);
            else
                pos->mergeTimers(it);
        }
    }

    timer() : name("no-name") {
    }

    void init() {
        isInitialized = true;
        prevTimer = activeTimer();
        activeTimer() = this;

        if (prevTimer == nullptr) {
            return;
        }

        std::vector<timer>::iterator begin = prevTimer->lowerTimers.begin();
        std::vector<timer>::iterator end = prevTimer->lowerTimers.end();
        std::vector<timer>::iterator pos = std::find(begin, end, *this);
        // duration = std::chrono::duration<double>(this->end - this->start).count();

        if (pos == end) {
            timer tmp = *this;
            tmp.isFinished = true;
            prevTimer->lowerTimers.push_back(tmp);
        }
    }

    static std::chrono::high_resolution_clock::time_point now() {
        return std::chrono::high_resolution_clock::now();
    }

    const char* name;
    bool isFinished{false};
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    double duration{0.0};
    int64_t calls{0};
    bool isInitialized{false};

    std::vector<timer> lowerTimers;

    timer* prevTimer{nullptr};

    static timer*& activeTimer() {
        static timer* value{nullptr};
        return value;
    }
};

extern timer globalTimer;

static inline void printTreeProfile(std::ostream& os = std::cout) {
    globalTimer.finish();
    globalTimer.printTimers(0, os);
}

void printTableFromTreeImpl(bool isHead, std::ostream& os, const std::vector<std::string_view>& names);
void printSliceImpl(std::ostream& os, const std::vector<std::string_view>& names);

template <typename... Types>
static inline void printTableFromTree(bool isHead, std::ostream& os, Types... names) {
    // if constexpr (sizeof...(Types) != 0)
    printTableFromTreeImpl(isHead, os, {names...});
}

template <typename... Types>
static inline void printSlice(std::ostream& os, Types... names) {
    printSliceImpl(os, {names...});
}

static inline void clearTree() {
    globalTimer.clear();
}

/*
 * Flat timer
 */

struct Event {
    const char* name;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    int64_t m;
    MeasureUnit unit;
    bool isOmp;
};

struct alignas(64) AlignedInt {
    int64_t val;
};

constexpr int64_t maxEvents = 1024 * 1024 * 1024 / sizeof(Event);
constexpr int64_t maxThreads = 128;
constexpr int64_t maxEventsPerThread = maxEvents / maxThreads;

extern Event events[maxEvents];
extern AlignedInt currEvents[maxThreads];

// global zero point for flat timers
extern std::chrono::high_resolution_clock::time_point globalStart;

class flatTimer {
   public:
    flatTimer(const char* nameIn, int64_t mIn = -1, MeasureUnit unitIn = MeasureUnit::not_set) {
        start(nameIn, mIn, unitIn);
    }

    // created for handling destructors calling
    flatTimer(NoStart) {
    }

    ~flatTimer() {
        finish();
    }

    void start(const char* nameIn, int64_t mIn = -1, MeasureUnit unitIn = MeasureUnit::not_set) {
        const int64_t thrnum = omp_get_thread_num();
        if (thrnum >= maxThreads) {
            isActive = false;
            return;
        }

        const int64_t currNum = currEvents[thrnum].val;
        if (currNum + 1 >= maxEventsPerThread) {
            name = "record limit per thread";
            eventNumber = maxEventsPerThread * thrnum + maxEventsPerThread - 1;
            static_assert(maxEventsPerThread > 1);
            startTime = events[eventNumber - 1].end;
            currEvents[thrnum].val = maxEventsPerThread;
        } else {
            eventNumber = maxEventsPerThread * thrnum + currNum;
            name = nameIn;
            m = mIn;
            unit = unitIn;

            currEvents[thrnum].val += 1;
            startTime = now();
        }
    }

    void finish() {
        if (!isActive) {
            return;
        }
        events[eventNumber].end = now();
        isActive = false;
        events[eventNumber].name = name;
        events[eventNumber].start = startTime;
        events[eventNumber].m = m;
        events[eventNumber].unit = unit;
        events[eventNumber].isOmp = omp_in_parallel();
    }

   private:
    static std::chrono::high_resolution_clock::time_point now() {
        return std::chrono::high_resolution_clock::now();
    }

    const char* name{};
    std::chrono::high_resolution_clock::time_point startTime;
    int64_t m{-1};
    MeasureUnit unit;
    int64_t eventNumber{};
    bool isActive = true;
};

class commonTimer {
   public:
    commonTimer(const char* nameIn, int64_t mIn = -1, MeasureUnit unitIn = MeasureUnit::not_set)
        : flat(nameIn, mIn, unitIn), tree(nameIn) {
    }

    void finish() {
        flat.finish();
        tree.finish();
    }

   private:
    flatTimer flat;
    timer tree;
};

extern void writeFullProfile(const char* filename);
extern void writeFullProfile(const std::string&);
extern void clearFullProfile();
}   // namespace timer

#define RECORD_TIMER_PARAMS(...)                                          \
    timer::timer _timer(std::source_location::current().function_name()); \
    timer::flatTimer _flatTimer(std::source_location::current().function_name(), __VA_ARGS__)

#define RECORD_TIMER                                                      \
    timer::timer _timer(std::source_location::current().function_name()); \
    timer::flatTimer _flatTimer(std::source_location::current().function_name())

#define RECORD_TIMER_NAMED(NAME) \
    timer::timer _timer(NAME);   \
    timer::flatTimer _flatTimer(NAME)

#else   // USE_TIMERS

namespace timer {

class timer {
   public:
    timer(const char*) {
    }

    void finish() {
    }

    ~timer() {
    }
};

class flatTimer {
   public:
    flatTimer(const char* nameIn, int64_t mIn = -1, MeasureUnit unitIn = MeasureUnit::not_set) {
        (void) nameIn;
        (void) mIn;
        (void) unitIn;
    }

    // created for handling destructors calling
    flatTimer(NoStart) {
    }

    ~flatTimer() {
        finish();
    }

    void start(const char* nameIn, int64_t mIn = -1, MeasureUnit unitIn = MeasureUnit::not_set) {
        (void) nameIn;
        (void) mIn;
        (void) unitIn;
    }

    void finish() {
    }
};

class commonTimer {
   public:
    commonTimer(const char* nameIn, int64_t mIn = -1, MeasureUnit unitIn = MeasureUnit::not_set) {
        (void) nameIn;
        (void) mIn;
        (void) unitIn;
    }
    void finish() {
    }
};

static inline void writeFullProfile(const char*) {};
static inline void writeFullProfile(const std::string&) {};
static inline void clearFullProfile() {};
static inline void printTreeProfile(std::ostream& os = std::cout) {
    (void) os;
};

static inline void clearTree() {
}

static inline void clear() {
}

template <typename... Types>
static inline void printTableFromTree(bool isHead, std::ostream& os, Types... names) {
    (void) isHead;
    (void) os;
    ((void) names, ...);
}

template <typename... Types>
static inline void printSlice(std::ostream& os, Types... names) {
    (void) os;
    ((void) names, ...);
}

#define RECORD_TIMER_PARAMS(...)
#define RECORD_TIMER
#define RECORD_TIMER_NAMED(NAME)

}   // namespace timer

#endif   // USE_TIMERS
