#pragma once

#include <cstdint>
#include <random>

/* While C++ random generators has linear complexity of discard function, this has logarithmic complexity: allows to use
 * it for reproducibility without perf degradation. Inspired by
 * https://lemire.me/blog/2019/03/19/the-fastest-conventional-random-number-generator-that-can-pass-big-crush/
 * https://github.com/lemire/testingRNG/blob/master/source/lehmer64.h
 */
struct LehmerEngine {
    using result_type = uint64_t;

    LehmerEngine() : state(1) {
    }

    LehmerEngine(uint64_t seed) {
        state = seed;   // zero state is degenerate for Lehmer64 (always returns 0)
        state += 1;     // avoid overflow
        discard(seed);
        discard(2);   // avoid overflow too;
    }

    uint64_t operator()() {
        state *= magicConstant;
        return state >> 64;
    }

    template <typename T>
    inline void discard(T count) {
        static_assert(std::is_integral_v<T>);
        if (count < 0) {
            discard<__uint128_t>(std::numeric_limits<__uint128_t>::max() - static_cast<__uint128_t>(-count) + 1);
            return;
        }
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
