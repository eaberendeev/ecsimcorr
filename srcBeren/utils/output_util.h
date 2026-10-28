// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef OUTPUT_UTIL_H_
#define OUTPUT_UTIL_H_

#include <string>
#include <vector>

#include "containers.h"

void write_field_to_file(const std::string& dataName, const Field3d& field);
void read_field_from_file(const std::string& dataName, Field3d& field);

#endif
