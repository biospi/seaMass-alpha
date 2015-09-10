set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 "1.8.5" NO_MODULE COMPONENTS C CXX)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 "1.8.5" COMPONENTS C CXX REQUIRED)
endif()

find_package(Boost 1.54 REQUIRED COMPONENTS
	program_options
)

add_executable(smpeak
	smpeak.cpp
	peakcore.cpp
	SMPFile.cpp
)
target_include_directories(smpeak PUBLIC
	${HDF5_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
)
target_link_libraries(smpeak
	${HDF5_LIBRARIES}
	${Boost_LIBRARIES}
)