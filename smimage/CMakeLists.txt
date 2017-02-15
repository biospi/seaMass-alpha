set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 "1.8.5" NO_MODULE COMPONENTS C CXX)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 "1.8.5" COMPONENTS C CXX REQUIRED)
endif()

find_package(Boost 1.54 REQUIRED COMPONENTS
	program_options
)

set(HEADER_FILES
	core.hpp
	imagecore.hpp
	#SMGWriter.hpp
	SmrWriter.hpp
)

set(SOURCE_FILES
	core.cpp
	imagecore.cpp
	#SMGWriter.cpp
	smimage.cpp
	SmrWriter.cpp
)

add_executable(smimage ${HEADER_FILES} ${SOURCE_FILES}
)
message(WTF !!! ${HDF5_INCLUDE_DIRS})
message(WTF !!! ${HDF5_LIBRARIES})

target_include_directories(seamass PUBLIC
	${HDF5_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
)
target_link_libraries(smimage LINK_PUBLIC
	seamass_kernel
	${HDF5_LIBRARIES}
	${Boost_LIBRARIES}
)
