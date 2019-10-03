cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(IPCpostprocess)

set(HELPERS
    ../helpers/randomNumberGenerator.hpp
    ../helpers/isotropic_pair_correlation_function.hpp
)
set(POSTPROCESS
    IPCpostprocess.cpp
    IPCpostprocess.hpp
    main.cpp
)
set(ALL_FILES
    ${HELPERS}
     ${POSTPROCESS}
    ../IPC.hpp
)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")

add_executable(${PROJECT_NAME} ${ALL_FILES})
