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

#for r better performance of function Mesh::update_Lmat2
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(mesh/Mesh.cpp PROPERTIES COMPILE_FLAGS "-Ofast -funroll-loops")
endif()

add_library(Fields OBJECT ${SOURCE_LIB})
