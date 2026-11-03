#pragma once

#include <nlohmann/json.hpp>

#include "fields/external_fieldsE.h"

// Фабрика для ExternalFieldE: принимает либо объект, либо массив объектов
// Поддерживается:
//  - { "uniform_field": { "value": [ex,ey,ez] } }
//  - { "uniformly_charged_cylinder": { "radius": R, "value": V } }
//  - или массив [{...}, {...}]
std::unique_ptr<ElectricFieldConfig> create_electric_field_config(const nlohmann::json& system_config,
                                                                  const char* key = "ExternalFieldE");
