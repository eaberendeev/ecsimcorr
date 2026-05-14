#ifndef DAMPING_H_
#define DAMPING_H_
#include "World.h"
#include "containers.h"
#include "nlohmann/json.hpp"
using json = nlohmann::json;
void Damping_Func(double& source, double i, double maxi, double& energyDamp);
double damping_fields_circleXY(Field3d& fieldE, Field3d& fieldB, const Domain& domain, const json& config);
double damping_fields_rectangle(Field3d& fieldE, Field3d& fieldB, const json& config);
double damping_fields(Field3d& fieldE, Field3d& fieldB, const Domain& domain, const json& config);
#endif
