find_package(netCDF NO_MODULE)
if(NOT netCDF_FOUND)
    find_package(OldNetCDF REQUIRED)
endif()

add_subdirectory(intel)

add_library(seamass_kernel
        kernel.hpp
        VecMat.hpp
        VecMat.tpp
        FileNetcdf.cpp
        FileNetcdf.hpp
        FileNetcdf.tpp
        )
target_include_directories(seamass_kernel PUBLIC
        ${netCDF_INCLUDE_DIR}
        )
target_link_libraries(seamass_kernel
        seamass_kernel_intel # could replace with nvidia implementation
        ${netCDF_LIBRARIES}
        )
