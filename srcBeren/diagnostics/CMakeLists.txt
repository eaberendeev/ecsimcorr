project(Diagnostics LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  DiagnosticEnergy.cpp
                Diagnostic.cpp 
                Timer.cpp)

add_library(Diagnostics OBJECT ${SOURCE_LIB})
