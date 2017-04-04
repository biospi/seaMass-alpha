find_package(Boost 1.54 COMPONENTS filesystem REQUIRED)
find_package(PugiXML REQUIRED)

set(HEADER_FILES
	Dataset.cpp
	DatasetMzmlb.hpp
	mzMLxml.hpp
	MzMLb.hpp
	NetcdfWriter.hpp
)
set(SOURCE_FILES
	Dataset.cpp
	DatasetMzmlb.cpp
	mzMLxml.cpp
	MzMLb.cpp
	NetcdfWriter.cpp
)
set(TEMPLATE_FILES
	mzMLxml.tpp
	DatasetMzmlb.tpp
)

source_group("Template Files" FILES ${TEMPLATE_FILES})
add_library(seamass_io ${HEADER_FILES} ${SOURCE_FILES} ${TEMPLATE_FILES})
target_include_directories(seamass_io PUBLIC
	${Boost_INCLUDE_DIR}
	${PugiXML_INCLUDE_DIR}
)
target_link_libraries(seamass_io
	seamass_core
	seamass_kernel
	${Boost_LIBRARY}
	${PugiXML_LIBRARIES}
)

