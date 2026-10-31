project(Fields LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  maxwell_operators.cpp
                faraday.cpp
                Coil.cpp
                Damping.cpp
                external_fieldsE.cpp
                external_fieldsB.cpp
                solverSLE.cpp)

# Faraday's update uses the same arithmetic optimization as the field kernels.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(faraday.cpp PROPERTIES COMPILE_FLAGS "-Ofast -funroll-loops")
endif()

add_library(Fields OBJECT ${SOURCE_LIB})
