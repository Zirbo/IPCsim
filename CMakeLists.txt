cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(IPCsim)

set(HELPERS_FILES
    src/cell_lists.hpp
    src/cell_lists.cpp
    src/randomNumberGenerator.hpp
)
set(IPC_CLASS_FILES
    src/IPC.hpp
    src/IPCsimulation.cpp
    src/IPCsimulation.hpp
)
set(MAIN_FILES
    src/main.cpp
)
set(ALL_FILES
    ${HELPERS_FILES}
    ${IPC_CLASS_FILES}
    ${MAIN_FILES}
)

add_executable(${PROJECT_NAME} ${ALL_FILES})
#target_compile_options(${PROJECT_NAME} PRIVATE "-O2 -g")
target_compile_options(${PROJECT_NAME} PRIVATE "-fexceptions")
set(CMAKE_BUILD_TYPE Debug)
