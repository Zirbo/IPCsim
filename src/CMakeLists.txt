cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(IPCsim)

set(HELPERS
    helpers/cell_lists.hpp
    helpers/cell_lists.cpp
    helpers/randomNumberGenerator.hpp
    helpers/isotropic_pair_correlation_function.hpp
)
set(SIMULATION_FILES
    IPC.hpp
    simulation/main.cpp
    simulation/IPCsimulation.hpp
    simulation/IPCsimulationCommon.cpp
    simulation/IPCsimulation.cpp
    simulation/IPCsimulationJanus.cpp
)
set(ALL_FILES
    ${HELPERS}
    ${SIMULATION_FILES}
)

find_package(OpenMP REQUIRED)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")

add_executable(${PROJECT_NAME} ${ALL_FILES})
target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_LIBRARIES})
configure_file(simulation/input.in input.in)
configure_file(simulation/stages.in stages.in)

message(STATUS "COMMON FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "DEBUG FLAGS: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "RELEASE FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
