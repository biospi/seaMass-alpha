find_package(Boost 1.54 COMPONENTS filesystem REQUIRED)
find_package(PugiXML REQUIRED)
find_package(SpatialIndex 1.8.5 REQUIRED)

set(HEADER_FILES
	MSFileData.hpp
	mzMLxml.hpp
	RTreeReader.hpp
	MzMLb.hpp
	NetcdfWriter.hpp
)
set(SOURCE_FILES
	MSFileData.cpp
	mzMLxml.cpp
	RTreeReader.cpp
	MzMLb.cpp
	NetcdfWriter.cpp
)
set(TEMPLATE_FILES
	mzMLxml.tpp
)

source_group("Template Files" FILES ${TEMPLATE_FILES})
add_library(seamass_io ${HEADER_FILES} ${SOURCE_FILES} ${TEMPLATE_FILES})
target_include_directories(seamass_io PUBLIC
	${PugiXML_INCLUDE_DIR}
	${SPATIALINDEX_INCLUDE_DIR}
	${Boost_INCLUDE_DIR}
)
target_link_libraries(seamass_io
	seamass_core
	seamass_kernel
	${PugiXML_LIBRARIES}
	${SPATIALINDEX_LIBRARY}
	${Boost_LIBRARY}
)

