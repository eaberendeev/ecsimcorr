project(Fields LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  mesh/operators.cpp
                mesh/Mesh.cpp 
                Coil.cpp
                Damping.cpp
                external_fieldsE.cpp
                external_fieldsB.cpp
                solverSLE.cpp
                )

add_library(Fields OBJECT ${SOURCE_LIB})
