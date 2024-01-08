find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")

find_package(IPP CONFIG REQUIRED PATHS $ENV{IPPROOT})
if(IPP_FOUND)
message(STATUS "oneIPP found")
endif()

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
        $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:IPP::ippcore,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:IPP::ipps,INTERFACE_INCLUDE_DIRECTORIES>
        )
target_link_libraries(seamass_kernel PUBLIC
        $<LINK_ONLY:MKL::MKL>
        $<LINK_ONLY:IPP::ippcore>
        $<LINK_ONLY:IPP::ipps>
        )
