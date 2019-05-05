set (CMAKE_CXX_STANDARD 11)

project(IPCsim)

set(HELPERS_FILES
  src/cell_lists.hpp
  src/cell_lists.cpp
  src/zilrandom.hpp
  src/zilvectors.hpp
)
set(IPC_CLASS_FILES
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
