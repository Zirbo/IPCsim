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


#find_package(OpenMP REQUIRED)
find_package(OpenACC REQUIRED)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenACC_CXX_FLAGS}")

add_executable(${PROJECT_NAME} ${ALL_FILES})
target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_LIBRARIES})

#message(STATUS "COMMON FLAGS: ${CMAKE_CXX_FLAGS}")
#message(STATUS "DEBUG FLAGS: ${CMAKE_CXX_FLAGS_DEBUG}")
#message(STATUS "RELEASE FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")
