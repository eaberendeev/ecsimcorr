// collisions_with_neutrals.cpp

#include "collision_processing.hpp"
#include "collision_handler.hpp"
#include "cross_section.hpp"
#include "utils.hpp"
#include "Vec.h"
#include <tuple>

// Функция моделирования столкновения заряженной частицы с нейтральной
// по алгоритму Vahedi and Surendra с нулевыми столкновениями

// На вход принимает: vcp, vn скорости зар. частицы и нейтрала в с mc

std::tuple<bool, double3, double3> collision_with_neutral(
    double3& vcp, double3& vn,
    double mcp, double mn,
    double ncp, double nn,
    double freq_max // Теперь передается извне
) {
    // Определяем тип заряженной частицы
    bool is_electron = (mcp == 1.0);
    bool is_proton = (mcp == 1836.0);

    if (!is_electron && !is_proton) {
        throw std::invalid_argument("Неизвестный тип заряженной частицы");
    }

    // Вычисляем относительную скорость
    double3 v_rel = vcp - vn;
    double E = compute_energy(v_rel, mcp); // Энергия в mc^2 (нормирована)
    
    // Вычисляем вероятность столкновения с учетом нулевых столкновений
    double P_null_collision = compute_collision_probability(freq_max);
    std::cout << P_null_collision;
    
    // Проверяем, произошло ли столкновение
    if (!check_collision(P_null_collision)) {
       return {false, vcp, vn};
    }
  
    // Выбираем тип столкновения на основе относительных вероятностей
    CollisionType collision_type = select_collision_type(E, mcp, nn, freq_max);
    
    // Обрабатываем столкновение и пересчитываем скорости
    return process_collision(collision_type, vcp, vn, mcp, mn);
}

int main2(){
    std::cout << "Hello, World!" << std::endl;
    double3 v(0.1, 0.1, 0.1); 
    collision_with_neutral(
    v, v,
     1., 1836.,
     1,  1,
    0.01 // Теперь передается извне
);
    
    return 0;
}