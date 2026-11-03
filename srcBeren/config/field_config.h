#pragma once

#include <nlohmann/json.hpp>

#include "fields/Damping.h"
#include "fields/external_fieldsB.h"
#include "fields/external_fieldsE.h"

// Фабрика для ExternalFieldE: принимает либо объект, либо массив объектов
// Поддерживается:
//  - { "uniform_field": { "value": [ex,ey,ez] } }
//  - { "uniformly_charged_cylinder": { "radius": R, "value": V } }
//  - или массив [{...}, {...}]
std::unique_ptr<ElectricFieldConfig> create_electric_field_config(const nlohmann::json& system_config,
                                                                  const char* key = "ExternalFieldE");

// Фабрика для ExternalFieldB
// Поддерживается:
//  - { "uniform_field": { "value": [bx,by,bz] } }
//  - { "coils": [ {z,R,I}, ... ] }
//  - массив [{...}, {...}] для явного порядка
std::unique_ptr<MagneticFieldConfig> create_magnetic_field_config(const nlohmann::json& system_config,
                                                                  const char* key = "ExternalFieldB");

// Read once before stepping; restart reconstructs the same layers from the input.
DampingParameters damping_parameters_from_config(const nlohmann::json& config);
