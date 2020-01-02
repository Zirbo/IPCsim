cmake_minimum_required(VERSION 3.5)

project(IPCsim)

set(ALL_FILES
    ../IPC.hpp
    main.cpp
    IPCsimulation.hpp
    IPCsimulationCommon.cpp
    IPCsimulation.cpp
    IPCsimulationJanus.cpp
    IPCsimulationPrint.cpp
    IPCsimulationAnalysis.cpp
)

add_executable(${PROJECT_NAME} ${ALL_FILES})
target_link_libraries(${PROJECT_NAME} helpers ${OpenMP_CXX_LIBRARIES})
configure_file(input.in input.in)
configure_file(stages.in stages.in)
