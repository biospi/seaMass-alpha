set(USE_MKL_64BIT_LIB On)
find_package(Intel REQUIRED)

set(HEADER_FILES
    Matrix.hpp
	MatrixMKL.hpp
	MatrixSparseMKL.hpp
)
set(SOURCE_FILES
	MatrixMKL.cpp
	MatrixSparseMKL.cpp
)

add_library(seamass_math ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(seamass_math PUBLIC
  ${MKL_INCLUDE_DIR}
  ${IPP_INCLUDE_DIR}
)
target_link_libraries(seamass_math
  ${INTEL_LIBRARIES}
)
