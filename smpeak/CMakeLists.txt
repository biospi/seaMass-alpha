set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 "1.8.5" NO_MODULE COMPONENTS C CXX)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 "1.8.5" COMPONENTS C CXX REQUIRED)
endif()

find_package(Boost 1.54 REQUIRED COMPONENTS
	program_options
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(PugiXML REQUIRED)

set(HEADER_FILES
	SMPFile.hpp
	BsplineData.hpp
	MathOperator.hpp
	PeakData.hpp
	PeakManager.hpp
	PeakOperator.hpp
	SMData.hpp
)
set(SOURCE_FILES
	smpeak.cpp
	SMPFile.cpp
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
