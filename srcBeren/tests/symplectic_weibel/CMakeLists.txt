# Explicitly invoked by tests/integration/Weibel_symplectic/run.py. This
# quantitative physical validation is too expensive for the unit suite.
add_executable(symplectic_weibel main.cpp)
target_link_libraries(symplectic_weibel ${LIBRARIES})
