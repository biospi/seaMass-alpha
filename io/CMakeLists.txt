set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 1.8.5 NO_MODULE)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 1.8.5 REQUIRED)
endif()

find_package(netCDF NO_MODULE)
if(NOT netCDF_FOUND)
    find_package(OldNetCDF REQUIRED)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(PugiXML REQUIRED)
find_package(SpatialIndex 1.8.5 REQUIRED)

set(HEADER_FILES
	MSFileData.hpp
    NetCDFile.hpp
	peakcore.hpp
)
set(SOURCE_FILES
	MSFileData.cpp
    NetCDFile.cpp
)
set(TEMPLATE_FILES
	NetCDFile.tpp
	peakcore.tpp
)

source_group("Template Files" FILES ${TEMPLATE_FILES})
add_library(seamass_io ${HEADER_FILES} ${SOURCE_FILES} ${TEMPLATE_FILES})
target_include_directories(seamass_io PUBLIC
	${HDF5_INCLUDE_DIR}
	${netCDF_INCLUDE_DIR}
	${PugiXML_INCLUDE_DIR}
	${SPATIALINDEX_INCLUDE_DIR}
)
target_link_libraries(seamass_io
	seamass_core
	${HDF5_LIBRARIES}
	${netCDF_LIBRARIES}
	${PugiXML_LIBRARIES}
	${SPATIALINDEX_LIBRARY}
)

