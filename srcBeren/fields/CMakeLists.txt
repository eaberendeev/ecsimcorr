project(Fields LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  solverSLE.cpp
                operators.cpp 
                Mesh.cpp 
                Damping.cpp
                Coil.cpp)

add_library(Fields OBJECT ${SOURCE_LIB})
