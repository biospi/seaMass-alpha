find_package(netCDF NO_MODULE)
if(NOT netCDF_FOUND)
    find_package(OldNetCDF REQUIRED)
endif()

set(USE_MKL_64BIT_LIB On)
find_package(Intel REQUIRED)

set(HEADER_FILES
    Matrix.hpp
	MatrixMKL.hpp
	MatrixSparseMKL.hpp
    VecMat.hpp
    NetcdfFile.hpp
)
set(SOURCE_FILES
	MatrixMKL.cpp
	MatrixSparseMKL.cpp
    NetcdfFile.cpp
)
set(TEMPLATE_FILES
    mzMLxml.tpp
    VecMat.tpp
    NetcdfFile.tpp
)

add_library(seamass_kernel ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(seamass_kernel PUBLIC
    ${netCDF_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR}
    ${IPP_INCLUDE_DIR}
)
target_link_libraries(seamass_kernel
    ${netCDF_LIBRARIES}
    ${INTEL_LIBRARIES}
)
