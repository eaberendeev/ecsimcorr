// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com

#pragma once

#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <omp.h>

#include <cmath>
#include <random>

class RandomGenerator {
   public:
    double Uniform01() { return urd(generator); };

    void SetRandSeed(int val) { generator.seed(val); }

    double Gauss(double sigma) {
        double r1 = Uniform01();
        double r2 = Uniform01();

        return sigma * sqrt(-2.0 * log(r1)) * sin(2.0 * M_PI * r2);
    }
    std::mt19937& gen() { return generator; }
    std::mt19937 generator;
    std::uniform_real_distribution<> urd{0, 1};
};

class ThreadRandomGenerator {
   public:
    ThreadRandomGenerator() : generators(omp_get_max_threads()) {
        for (size_t i = 0; i < generators.size(); i++) {
            generators[i].SetRandSeed(100 * i + 1);
        }
    }
    double Uniform01() {
        int i = omp_get_thread_num();
        return generators[i].Uniform01();
    };
    double Gauss(double sigma) {
        int i = omp_get_thread_num();
        return generators[i].Gauss(sigma);
    };
    void SetRandSeed(int val) {
        int i = omp_get_thread_num();
        generators[i].SetRandSeed(val+i);
    };
    std::mt19937& gen() {
        int i = omp_get_thread_num();
        return generators[i].gen();
    }

   private:
    std::vector<RandomGenerator> generators;
};

#endif   // RANDOM_GENERATOR_H
