cmake_minimum_required(VERSION 3.7)
project(tomography)

set(CMAKE_CXX_STANDARD 17)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

set(SOURCE_FILES
        src/dimension.cpp
        src/distribution.cpp
        src/grid.cpp
        src/integration.cpp
        src/limits.cpp
        src/main.cpp
        src/point.cpp
        src/ray.cpp
        src/reconstruction.cpp
        src/solution.cpp
        src/spline.cpp
        src/tools.cpp
        src/transformation.cpp
        src/vector_sparse.cpp
        )
add_executable(tomography ${SOURCE_FILES})

target_link_libraries(tomography "-lstdc++fs")
