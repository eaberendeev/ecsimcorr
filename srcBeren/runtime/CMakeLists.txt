add_library(beren_runtime OBJECT profiler.cpp)
target_link_libraries(beren_runtime PUBLIC beren_core OpenMP::OpenMP_CXX)

# Only executables linking the wrapper implementation request linker interposition.
# Pure geometry/configuration clients need neither wrappers nor runtime objects.
if(MEMORY_TIMERS_LEVEL GREATER 0)
    message(STATUS "Common memory profiling")
    target_link_options(beren_runtime INTERFACE
        -Wl,--wrap=realloc -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=reallocarray
        -Wl,--wrap=free -Wl,--wrap=_Znwm -Wl,--wrap=_Znam -Wl,--wrap=_ZdlPvm)
endif()
if(MEMORY_TIMERS_LEVEL GREATER 1)
    message(STATUS "Extended memory profiling")
    target_link_options(beren_runtime INTERFACE -Wl,--wrap=memcpy -Wl,--wrap=memset)
endif()
