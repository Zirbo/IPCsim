cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

project(IPCsim)

set(HELPERS
    src/cell_lists.hpp
    src/cell_lists.cpp
    src/randomNumberGenerator.hpp
)
set(COMMON_FILES
    ${HELPERS}
    src/IPC.hpp
)
set(IPC
    src/IPCsimulation.cpp
    src/IPCsimulation.hpp
    src/main.cpp
)
set(JANUS
    src/JanusIPCsimulation.cpp
    src/JanusIPCsimulation.hpp
    src/mainJanus.cpp
)

find_package(OpenMP REQUIRED)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} bld)

add_executable(IPCsim ${COMMON_FILES} ${IPC})
target_link_libraries(IPCsim ${OpenMP_CXX_LIBRARIES})

add_executable(JanusIPCsim ${COMMON_FILES} ${JANUS})
target_link_libraries(JanusIPCsim ${OpenMP_CXX_LIBRARIES})

message(STATUS "COMMON FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "DEBUG FLAGS: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "RELEASE FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
