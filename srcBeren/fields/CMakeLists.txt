project(fields LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  maxwell_operators.cpp
                field_boundaries.cpp
                Coil.cpp
                Damping.cpp
                external_fieldsE.cpp
                external_fieldsB.cpp
                solverSLE.cpp)

add_library(fields OBJECT ${SOURCE_LIB})
