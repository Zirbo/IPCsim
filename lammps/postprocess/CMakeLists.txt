cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")


project(lammpsIPCpostprocess)

set(ALL_FILES
    IPC.hpp
    IPCpostprocess.cpp
    IPCpostprocess.hpp
    main.cpp
    IPCpostprocessPotential.hpp
    IPCpostprocessPotential.cpp
   # IPCpostprocessDataAnalysis.hpp
   # IPCpostprocessDataAnalysis.cpp
)

add_executable(${PROJECT_NAME} ${ALL_FILES})
