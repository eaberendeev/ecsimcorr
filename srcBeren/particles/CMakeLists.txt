project(Particles LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  collision.cpp
                Particles.cpp 
                ParticlesDistribution.cpp 
                ParticlesDiagnostic.cpp 
                ParticlesMethods.cpp)

add_library(Particles OBJECT ${SOURCE_LIB})
