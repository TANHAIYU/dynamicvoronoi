project(DYNAMICVORONOI)
cmake_minimum_required(VERSION 3.0)

SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Release)
find_package (Eigen3 REQUIRED )

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib )
SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib )


# COMPILER FLAGS
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")

include_directories(${EIGEN3_INCLUDE_DIR})
add_library(dynamicvoronoi src/dynamicvoronoi.cpp)

add_executable(example src/example.cpp)
target_link_libraries(example dynamicvoronoi Eigen3::Eigen)

