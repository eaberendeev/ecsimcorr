#pragma once

#include <nlohmann/json.hpp>

#include "fields/external_fieldsB.h"

// Фабрика для ExternalFieldB
// Поддерживается:
//  - { "uniform_field": { "value": [bx,by,bz] } }
//  - { "coils": [ {z,R,I}, ... ] }
//  - массив [{...}, {...}] для явного порядка
std::unique_ptr<MagneticFieldConfig> create_magnetic_field_config(const nlohmann::json& system_config,
                                                                  const char* key = "ExternalFieldB");
