#pragma once

#include <tuple>

#include "Vec.h"
#include "cross_section.hpp"
#include "sgs.h"

class ColliderWithNeutrals {
   public:
    ColliderWithNeutrals(const double n0) {
        const double wp = SGS::get_plasma_freq(n0);
        const double l0 = SGS::c / wp;
        dpd0 = n0 / l0;
    }   // Функция моделирования столкновения заряженной частицы с нейтральной
        // по алгоритму Vahedi and Surendra с нулевыми столкновениями

    // На вход принимает: vcp, vn скорости зар. частицы и нейтрала в с mc

    std::tuple<bool, double3, double3> collision_with_neutral(
        double3& vcp, double3& vn, double mcp, double mn, double ncp, double nn, double dt,
        double freq_max   // Теперь передается извне
    );

   private:
    double dpd0;   // density per distance, n0/l0

    // Функции для вычисления сечений столкновений
    double Sigma_e(double E);   // Ударная ионизация электронами
    double Sigma_p(double E);   // Ударная ионизация протонами
    double Sigma_cx(double E);   // Резонансная перезарядка

    // Функция проверки, произошло ли столкновение
    bool check_collision(double P_collision);

    // Функция выбора типа столкновения
    CollisionType select_collision_type(double E, double mcp, double nn,
                                        double freq_max);
};