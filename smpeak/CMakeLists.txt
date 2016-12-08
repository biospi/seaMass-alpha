find_package(Boost 1.54 COMPONENTS program_options REQUIRED)
find_package(PugiXML REQUIRED)

set(HEADER_FILES
	BsplineData.hpp
	MathOperator.hpp
	PeakData.hpp
	PeakManager.hpp
	PeakOperator.hpp
	SMData.hpp
)
set(SOURCE_FILES
	smpeak.cpp
)
set(TEMPLATE_FILES
	BsplineData.tpp
	MathOperator.tpp
	PeakData.tpp
	PeakOperator.tpp
	SMData.tpp
)

source_group("Template Files" FILES ${TEMPLATE_FILES})
add_executable(smpeak ${HEADER_FILES} ${SOURCE_FILES} ${TEMPLATE_FILES})
target_include_directories(smpeak PUBLIC
	${HDF5_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${PugiXML_INCLUDE_DIR}
)
target_link_libraries(smpeak
	seamass_io
	seamass_core
	${HDF5_LIBRARIES}
	${Boost_LIBRARIES}
	${PugiXML_LIBRARIES}
)
