cmake_minimum_required(VERSION 2.8.12)
project(seamass C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(OpenMP REQUIRED)
set(CXX_EXTRA_FLAGS "-std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${CXX_EXTRA_FLAGS}")
find_package(IntelOpenMP REQUIRED)

add_subdirectory(seamass)
add_subdirectory(restoration)
add_subdirectory(smpeak)
