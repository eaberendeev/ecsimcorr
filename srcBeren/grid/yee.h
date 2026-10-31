#pragma once

#include "core/vector3.h"

enum class FieldType { ELECTRIC, MAGNETIC, DENSITY, CURRENT };

// Смещения для узлов полей относительно центра ячейки
inline Vector3R electric_shift(int component) {
    // Ex сдвинута на +dx/2 по X, Ey на +dy/2 по Y, Ez на +dz/2 по Z
    Vector3R s{0, 0, 0};
    if (component == 0)
        s.x() = 0.5;
    else if (component == 1)
        s.y() = 0.5;
    else if (component == 2)
        s.z() = 0.5;
    return s;
}

inline Vector3R magnetic_shift(int component) {
    // Bx сдвинута на +dy/2 и +dz/2, By на +dx/2 и +dz/2, Bz на +dx/2 и +dy/2
    assert(component >= 0 && component < 3);
    Vector3R s{0.5, 0.5, 0.5};
    s[component] = 0;
    return s;
}

inline Vector3R density_shift() {
    return {0, 0, 0};   // плотность в центре ячейки
}

// J and E share Yee edges in the discrete J dot E power density.
inline Vector3R yee_shift(FieldType field, int component) {
    if (field == FieldType::ELECTRIC || field == FieldType::CURRENT)
        return electric_shift(component);
    if (field == FieldType::MAGNETIC)
        return magnetic_shift(component);
    return density_shift();
}
