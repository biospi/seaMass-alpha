set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 "1.8.5" NO_MODULE COMPONENTS C CXX)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 "1.8.5" COMPONENTS C CXX REQUIRED)
endif()

add_executable(restoration
	restoration.cpp
)
target_include_directories(seamass PUBLIC
	${HDF5_INCLUDE_DIRS}
)
target_link_libraries (restoration LINK_PUBLIC
	seamass
	${HDF5_LIBRARIES}
)
