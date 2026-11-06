# Explicit experiment; never part of the fast CTest suite.
add_executable(weibel_whitney main.cpp)
target_link_libraries(weibel_whitney ${LIBRARIES})
# Separate binary so an experiment can finish with its original executable
# while a validated cell-assembly implementation is built alongside it.
add_executable(weibel_whitney_cell main.cpp)
target_link_libraries(weibel_whitney_cell ${LIBRARIES})
