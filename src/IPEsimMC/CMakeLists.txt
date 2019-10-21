cmake_minimum_required(VERSION 3.5)

project(IPEsim)

set(ALL_FILES
#    ../IPC.hpp
    main.cpp
    IPEsimulation.hpp
    IPEsimulation.cpp
)

add_executable(${PROJECT_NAME} ${ALL_FILES})
target_link_libraries(${PROJECT_NAME}) # helpers})
configure_file(input.in input.in)
configure_file(stages.in stages.in)
