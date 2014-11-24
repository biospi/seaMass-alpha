set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(SpatialIndex REQUIRED)
find_package(Boost 1.54 REQUIRED COMPONENTS system filesystem)

add_executable(viz_client
	viz_client.cpp
)
target_include_directories(viz_client PUBLIC
  ${Boost_INCLUDE_DIRS}
  ${SPATIALINDEX_INCLUDE_DIR}
)
target_link_libraries (viz_client LINK_PUBLIC
  seamass
  ${Boost_LIBRARIES}
  ${SPATIALINDEX_LIBRARY}
)
