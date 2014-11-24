cmake_minimum_required(VERSION 2.8.12)
project(seamass C CXX)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_subdirectory(seamass)
add_subdirectory(restoration)
add_subdirectory(viz_client)
