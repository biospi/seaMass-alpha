if(WIN32)
    find_package(HDF5 COMPONENTS C HL NO_MODULE REQUIRED static)
endif()

find_package(netCDF NO_MODULE)
if(NOT netCDF_FOUND)
    find_package(OldNetCDF REQUIRED)
endif()

add_library(seamass_io
        VecMat.hpp
        VecMat.tpp
        FileNetcdf.hpp
        FileNetcdf.cpp
        FileNetcdf.tpp
        ObserverMatrix.cpp
        ObserverMatrix.hpp
        ObserverMatrixSparse.cpp
        ObserverMatrixSparse.hpp
        )
target_include_directories(seamass_io PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${netCDF_INCLUDE_DIR}
        )
target_link_libraries(seamass_io
        seamass_kernel
        ${netCDF_LIBRARIES}
        )