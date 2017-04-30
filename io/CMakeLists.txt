find_package(netCDF NO_MODULE)
if(NOT netCDF_FOUND)
    find_package(OldNetCDF REQUIRED)
endif()

add_library(seamass_io
        FileNetcdf.hpp
        FileNetcdf.cpp
        FileNetcdf.tpp
        )
target_include_directories(seamass_io PUBLIC
        ${netCDF_INCLUDE_DIR}
        )
target_link_libraries(seamass_io
        seamass_kernel_intel # could replace with nvidia implementation
        ${netCDF_LIBRARIES}
        )