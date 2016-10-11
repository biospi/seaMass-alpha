find_package(Boost 1.54 COMPONENTS filesystem REQUIRED)
find_package(HDF5 COMPONENTS C CXX NO_MODULE REQUIRED static)
find_package(netCDF 4.4.1 REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(PugiXML REQUIRED)
find_package(SpatialIndex 1.8.5 REQUIRED)

set(HEADER_FILES
	MSFileData.hpp
    NetCDFile.hpp
	HDF5Writer.hpp
	iomath.hpp
	mzMLxml.hpp
	RTreeReader.hpp
	MzMLb.hpp
)
set(SOURCE_FILES
	MSFileData.cpp
    NetCDFile.cpp
	HDF5Writer.cpp
	mzMLxml.cpp
	RTreeReader.cpp
	MzMLb.cpp
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
	${Boost_INCLUDE_DIR}
)
target_link_libraries(seamass_io
	seamass_core
	${HDF5_LIBRARIES}
	${netCDF_LIBRARIES}
	${PugiXML_LIBRARIES}
	${SPATIALINDEX_LIBRARY}
	${Boost_LIBRARY}
)

