# To build examples provided in this book always use recommended commands:
#
# * cmake -B <build tree> -S <source tree>
# * cmake --build <build tree>
# 
# Specifying the minimum CMake version
cmake_minimum_required(VERSION 3.20)

## Defining languages and metadata
set(CURRENT_VERSION 1.0)

project(ra_solver
        VERSION ${CURRENT_VERSION}
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add local cmake module path
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# dependent packages
find_package(GUROBI REQUIRED)
find_package(cxxopts CONFIG REQUIRED)

# solving RAP with MILP models
add_executable(ra_mm
        src/mm_main.cpp
        src/mm_instance.cpp
        src/mm_instance.h
        src/mm_reliability.cpp
        src/mm_reliability.h
        src/mm_solution.cpp
        src/mm_solution.h
        src/mm_solver.cpp
        src/mm_solver.h)
target_include_directories(ra_mm PRIVATE ${GUROBI_INCLUDE_DIRS})
target_link_libraries(ra_mm optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
target_link_libraries(ra_mm ${GUROBI_LIBRARY} cxxopts::cxxopts)

# solving RAP with branch-and-bound
add_executable(ra_bb
        src/bb_main.cpp
        src/bb_instance.cpp
        src/bb_instance.h
        src/bb_solver.cpp
        src/bb_solver.h
        src/bb_utility.h
        src/bb_utility.cpp)
target_link_libraries(ra_bb PRIVATE cxxopts::cxxopts)



