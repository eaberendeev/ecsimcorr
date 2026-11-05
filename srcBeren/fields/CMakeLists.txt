add_library(fields OBJECT
    maxwell_operators.cpp
    field_boundaries.cpp
    Coil.cpp
    Damping.cpp
    external_fieldsE.cpp
    external_fieldsB.cpp
    solverSLE.cpp
)
target_link_libraries(fields PUBLIC beren_build_options beren_grid utils)
