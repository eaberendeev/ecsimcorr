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
        state = seed + 1;
        if (state == 0) {
            state = 1;   // zero state is degenerate for Lehmer64 (always returns 0)
        }
        discard(seed + 2);
    }

    uint64_t operator()() {
        state *= magicConstant;
        return state >> 64;
    }

    // advances the state by `count` steps in O(log count) via modular exponentiation;
    // unlike the previous bit-mask loop it cannot overflow for count >= 2^63
    void discard(int64_t count) {
        uint64_t c = static_cast<uint64_t>(count);
        __uint128_t mult = magicConstant;
        __uint128_t acc = 1;
        while (c != 0) {
            if (c & 1) {
                acc = acc * mult;
            }
            mult = mult * mult;
            c >>= 1;
        }
        state = state * acc;
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
