// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com

#pragma once

#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <random>
#include <cmath>

class RandomGenerator{
    public:
    std::mt19937 gen;
    std::uniform_real_distribution<> urd{0, 1};

    double Uniform01() { return urd(gen); };

    void SetRandSeed(int val) { gen.seed(val); }

    double Gauss(double sigma) {
        double r1 = Uniform01();
        double r2 = Uniform01();

        return sigma * sqrt(-2.0 * log(r1)) * sin(2.0 * M_PI * r2);
    }
};

#endif   // RANDOM_GENERATOR_H
