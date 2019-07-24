cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(IPCsim)

set(HELPERS
    ../helpers/cell_lists.hpp
    ../helpers/cell_lists.cpp
    ../helpers/randomNumberGenerator.hpp
)
set(COMMON_FILES
    ${HELPERS}
    ../IPC.hpp
)
set(IPC
    IPCsimulation.cpp
    IPCsimulation.hpp
    main.cpp
)

find_package(OpenMP REQUIRED)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")

add_executable(${PROJECT_NAME} ${COMMON_FILES} ${IPC})
target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_LIBRARIES})

message(STATUS "COMMON FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "DEBUG FLAGS: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "RELEASE FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
