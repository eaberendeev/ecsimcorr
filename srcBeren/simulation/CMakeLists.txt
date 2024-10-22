project(Simulation LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)

set(SOURCE_LIB  simulation.cpp
                recovery.cpp
                simulation_output.cpp)

add_library(Simulation OBJECT ${SOURCE_LIB})
