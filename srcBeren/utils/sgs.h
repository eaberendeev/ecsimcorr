// Author: Evgeny Berendeev

#pragma once

#ifndef SGS_H
#define SGS_H

#include <cmath>

namespace SGS {

constexpr double me = 9.10938356e-28;
constexpr double qe = 4.80320427e-10;
constexpr double c = 2.99792458e10;
constexpr double MC2 = 511.0;

inline double get_plasma_freq(double n0) {
    return pow(4 * M_PI * n0 * qe * qe / me, 0.5);
}

}   // namespace SGS

#endif
