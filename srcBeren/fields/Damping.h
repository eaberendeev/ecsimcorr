#ifndef DAMPING_H_
#define DAMPING_H_
#include "World.h"
#include "containers.h"
#include "nlohmann/json.hpp"
using json = nlohmann::json;

// Затухание одной компоненты поля в поглощающем слое.
// i     — расстояние от внешней кромки слоя (в ячейках),
// maxi  — толщина слоя,
// coeff — коэффициент затухания на внешней кромке (профиль квадратичный:
//         damp = coeff + (1 - coeff) * (i / maxi)^2, на внутренней границе -> 1).
void Damping_Func(double& source, double i, double maxi, double& energyDamp, double coeff);
double damping_fields_circleXY(Field3d& fieldE, Field3d& fieldB, const Domain& domain, const json& config,
                               double coeff);
double damping_fields_rectangle(Field3d& fieldE, Field3d& fieldB, const json& config, double coeff);

// Применяет поглощающие слои к полям; возвращает энергию полей, съеденную
// слоями. Тип слоя и толщины задаются в system_config: DampingType
// ("None"|"Rectangle"|"CircleXY"), DampCellsX/Y/Z_glob, DampingCoeff.
double damping_fields(Field3d& fieldE, Field3d& fieldB, const Domain& domain, const json& config);
#endif
