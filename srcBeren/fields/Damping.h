#ifndef DAMPING_H_
#define DAMPING_H_
#include <array>

#include "utils/containers.h"
#include "grid/domain.h"

enum class DampingType { None, Rectangle, CircleXY };

// Layer widths are in grid cells from the low/high outer edge of the stored mesh (including ghosts).
// CircleXY uses max(x_cells) and the X mesh spacing for its radial layer.
struct DampingParameters {
    DampingType type = DampingType::None;
    double coefficient = 0.8;
    std::array<int, 2> x_cells{};
    std::array<int, 2> y_cells{};
    std::array<int, 2> z_cells{};
};

// Затухание одной компоненты поля в поглощающем слое.
// i     — расстояние от внешней кромки слоя (в ячейках),
// maxi  — толщина слоя,
// coeff — коэффициент затухания на внешней кромке (профиль квадратичный:
//         damp = coeff + (1 - coeff) * (i / maxi)^2, на внутренней границе -> 1).
void Damping_Func(double& source, double i, double maxi, double& energyDamp, double coeff);
double damping_fields_circleXY(Field3d& fieldE, Field3d& fieldB, const Domain& domain,
                               const DampingParameters& parameters);
double damping_fields_rectangle(Field3d& fieldE, Field3d& fieldB, const DampingParameters& parameters);

// Применяет поглощающие слои к полям; возвращает энергию полей, съеденную
// слоями. Параметры постоянны в течение расчёта; JSON разбирается при инициализации.
double damping_fields(Field3d& fieldE, Field3d& fieldB, const Domain& domain, const DampingParameters& parameters);
#endif
