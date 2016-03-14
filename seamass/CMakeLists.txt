set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 1.8.5 NO_MODULE)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 1.8.5 REQUIRED)
endif()

find_package(MKL REQUIRED)
find_package(SpatialIndex 1.8.5 REQUIRED)
find_package(Boost 1.54 REQUIRED COMPONENTS
	system filesystem
)
find_package(netCDF REQUIRED)

add_library(seamass
  seamass.hpp
  OptimiserASRL.hpp
  SMOWriter.hpp
  SMVWriter.hpp
  BasisFunctions.hpp
  core.hpp
  seamass.cpp
  OptimiserASRL.cpp
  SMOWriter.cpp
  SMVWriter.cpp
  BasisFunctions.cpp
  core.cpp
  NetCDFile.hpp
  NetCDFile.tpp
  NetCDFile.cpp
)
target_include_directories(seamass PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MKL_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${SPATIALINDEX_INCLUDE_DIR}
  ${netCDF_INCLUDE_DIR}
  ../smpeak
)
target_link_libraries(seamass
  ${IntelOpenMP_LIBRARIES}
  ${MKL_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${Boost_LIBRARIES}
  ${SPATIALINDEX_LIBRARY}
  ${netCDF_LIBRARIES}
)

install(TARGETS seamass DESTINATION lib)
install(FILES seamass.hpp DESTINATION include)
