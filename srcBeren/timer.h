#pragma once

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>
#include <source_location>

namespace timer {

class timer {
   public:
    timer(std::string&& name) : name(name), start(now()) {
        if (!omp_in_parallel())
            init();
    }

    timer(const std::string& name) : name(name), start(now()) {
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
        return name == over.name;
    }

    void printTimers(const int64_t nestingDepth, std::ostream& os) const {
        for (int64_t i = 0; i < nestingDepth; ++i) os << "|  ";
        const size_t endPos = name.find('(');
        const std::string cuttedName = endPos == std::string::npos ? name : name.substr(0, endPos);
        os << "> " << cuttedName << ": " << calls << "[calls] " << duration << "[s]" << std::endl;
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

    void clear() {
        isFinished = false;
        duration = 0.0;
        calls = 1;
        lowerTimers.clear();
        start = now();
    }

   private:
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

    const std::string name;
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

inline void print(std::ostream& os = std::cout) {
    globalTimer.finish();
    globalTimer.printTimers(0, os);
}

inline void clear() {
    globalTimer.clear();
}

}   // namespace timer

#define RECORD_TIMER timer::timer _timer(std::source_location::current().function_name())