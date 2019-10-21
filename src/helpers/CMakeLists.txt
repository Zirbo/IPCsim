cmake_minimum_required(VERSION 3.5)

project(helpers)

set(ALL_FILES
    ../IPC.hpp
    cell_lists.hpp
    cell_lists.cpp
    randomNumberGenerator.hpp
    isotropic_pair_correlation_function.hpp
)

add_library(${PROJECT_NAME} ${ALL_FILES})
