cmake_minimum_required(VERSION 3.5)

project(IPCpostprocess)

set(ALL_FILES
    ../IPC.hpp
    IPCpostprocess.cpp
    IPCpostprocess.hpp
    main.cpp
)

add_executable(${PROJECT_NAME} ${ALL_FILES})
target_link_libraries(${PROJECT_NAME} helpers)
