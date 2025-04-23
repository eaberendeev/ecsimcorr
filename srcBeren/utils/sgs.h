// Author: Evgeny Berendeev

#pragma once

#ifndef SGS_H
#define SGS_H

#include <cmath>

namespace SGS {

const double me = 9.10938356e-28;
const double qe = 4.80320427e-10;
const double c = 2.99792458e10;
const double MC2 = 511.; 

inline double get_plasma_freq(double n0) {
    return pow(4 * M_PI * n0 * qe * qe / me, 0.5);
}

}   // namespace SGS

#endif
