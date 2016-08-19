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
	iomath.hpp
	mzMLxml.hpp
)
set(SOURCE_FILES
	MSFileData.cpp
    NetCDFile.cpp
	mzMLxml.cpp
)
set(TEMPLATE_FILES
	NetCDFile.tpp
	iomath.tpp
	mzMLxml.tpp
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

