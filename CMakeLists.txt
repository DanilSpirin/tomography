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
        tomography/dimension.cpp
        tomography/distribution.cpp
        tomography/grid.cpp
        tomography/integration.cpp
        tomography/main.cpp
        tomography/point.cpp
        tomography/ray.cpp
        tomography/reconctruction.cpp
        tomography/solution.cpp
        tomography/spline.cpp
        tomography/tools.cpp
        tomography/transformation.cpp
        tomography/vector_sparse.cpp
        )
add_executable(tomography ${SOURCE_FILES})

if (Boost_FOUND)
    target_link_libraries(tomography ${Boost_LIBRARIES})
endif()
