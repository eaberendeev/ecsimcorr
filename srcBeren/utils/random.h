#pragma once

#include <cstdint>
#include <random>

/* While C++ random generators has linear complexity of discard function, this has logarithmic complexity
 * Inspired by
 * https://lemire.me/blog/2019/03/19/the-fastest-conventional-random-number-generator-that-can-pass-big-crush/
 * https://github.com/lemire/testingRNG/blob/master/source/lehmer64.h
 */
struct LehmerEngine {
    using result_type = uint64_t;

    LehmerEngine() : state(1) {
    }

    LehmerEngine(uint64_t seed) {
        state = seed + 1;   // zero state is degenerate for Lehmer64 (always returns 0)
        discard(seed + 2);
    }

    uint64_t operator()() {
        state *= magicConstant;
        return state >> 64;
    }

    void discard(int64_t count) {
        __uint128_t mult = magicConstant;
        while (count != 0) {
            if (count & 1) {
                state = state * mult;
            }
            mult *= mult;
            count = count >> 1;
        }
    }

    static constexpr uint64_t min() {
        return 0;
    }

    static constexpr uint64_t max() {
        return ~(0UL);
    }

    __uint128_t state;

    static constexpr uint64_t magicConstant = 0xda942042e4dd58b5ull;
};
