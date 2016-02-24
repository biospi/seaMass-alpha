set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 "1.8.5" NO_MODULE COMPONENTS C CXX)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 "1.8.5" COMPONENTS C CXX REQUIRED)
endif()

find_package(Boost 1.54 REQUIRED COMPONENTS
	program_options
)

find_package(netCDF REQUIRED)

find_package(PugiXML REQUIRED)

add_executable(smpeak
	smpeak.cpp
	SMPFile.cpp
	NetCDFile.cpp
	BsplineData.hpp
	BsplineData.tpp
	MathOperator.hpp
	MathOperator.tpp
	NetCDFile.hpp
	NetCDFile.tpp
	peakcore.hpp
	peakcore.tpp
	PeakData.hpp
	PeakData.tpp
	PeakManager.hpp
	PeakOperator.hpp
	PeakOperator.tpp
	SMData.hpp
	SMData.tpp
	SMPFile.hpp
)
target_include_directories(smpeak PUBLIC
	${HDF5_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${netCDF_INCLUDE_DIR}
	${PugiXML_INCLUDE_DIR}
)
target_link_libraries(smpeak
	seamass
	${HDF5_LIBRARIES}
	${Boost_LIBRARIES}
	${netCDF_LIBRARIES}
	${PugiXML_LIBRARIES}
)
