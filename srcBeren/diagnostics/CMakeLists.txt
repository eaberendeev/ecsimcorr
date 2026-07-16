project(Diagnostics LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)

set(SOURCE_LIB  DiagnosticEnergy.cpp
                Diagnostic.cpp 
                Timer.cpp
                DiagnosticOutput.cpp)

add_library(Diagnostics OBJECT ${SOURCE_LIB})
