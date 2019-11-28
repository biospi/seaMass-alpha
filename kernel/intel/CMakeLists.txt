#set(USE_MKL_64BIT_LIB On)
find_package(Intel REQUIRED)

add_library(seamass_kernel
        ../Subject.cpp
        ../Subject.hpp
        ../SubjectMatrix.cpp
        ../SubjectMatrix.hpp
        ../SubjectMatrixSparse.cpp
        ../SubjectMatrixSparse.hpp
        ../Observer.cpp
        ../Observer.hpp
        types.hpp
        kernel.cpp
        kernel.hpp
        Matrix.cpp
        Matrix.hpp
        MatrixSparse.cpp
        MatrixSparse.hpp
        )
target_include_directories(seamass_kernel PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${MKL_INCLUDE_DIR}
        ${IPP_INCLUDE_DIR}
        )
target_link_libraries(seamass_kernel
        ${INTEL_LIBRARIES}
        )
