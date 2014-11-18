set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(MKL REQUIRED)

set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 "1.8.5" NO_MODULE)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 "1.8.5" REQUIRED)
endif()

add_library(seamass
  seamass.hpp
  OptimiserASRL.hpp
  HDF5File.hpp   
  BasisFunctions.hpp
  core.hpp
  seamass.cpp
  OptimiserASRL.cpp
  HDF5File.cpp   
  BasisFunctions.cpp
  core.cpp
)
target_include_directories(seamass PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MKL_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}
)
target_link_libraries(seamass
  ${MKL_LIBRARIES}
  ${HDF5_LIBRARIES}
)

install(TARGETS seamass DESTINATION lib)
install(FILES seamass.hpp DESTINATION include)
