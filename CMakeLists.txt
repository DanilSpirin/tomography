cmake_minimum_required(VERSION 3.7)
project(tomography)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost COMPONENTS system filesystem REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

if(Boost_FOUND)
    message("Found this boost libraries ${Boost_LIBRARIES}")
    message("Found this boost headers ${Boost_INCLUDE_DIR}")
    include_directories(${Boost_INCLUDE_DIR})
endif()

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

if (Boost_FOUND)
    target_link_libraries(tomography ${Boost_LIBRARIES})
endif()
