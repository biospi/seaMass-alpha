#set(USE_MKL_64BIT_LIB On)
find_package(Intel REQUIRED)

add_library(seamass_kernel_intel
        kernel.cpp
        kernel.hpp
        Matrix.cpp
        Matrix.hpp
        MatrixSparse.cpp
        MatrixSparse.hpp
        types.hpp
        Callback.cpp
        Callback.hpp
        )
target_include_directories(seamass_kernel_intel PUBLIC
        ${MKL_INCLUDE_DIR}
        ${IPP_INCLUDE_DIR}
        )
target_link_libraries(seamass_kernel_intel
        ${INTEL_LIBRARIES}
        )
