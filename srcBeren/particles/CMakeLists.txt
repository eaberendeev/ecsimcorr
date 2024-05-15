project(Particles LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  collision.cpp
                ParticlesArray.cpp 
                ParticlesDistribution.cpp 
                ParticlesDiagnostic.cpp 
                ParticlesCore.cpp)

add_library(Particles OBJECT ${SOURCE_LIB})
