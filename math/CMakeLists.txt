find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${CXX_EXTRA_FLAGS}")
endif ()

set(USE_MKL_64BIT_LIB On)
find_package(MKL REQUIRED)

set(HEADER_FILES
	Matrix.hpp
	MatrixSparse.hpp
)
set(SOURCE_FILES
	Matrix.cpp
	MatrixSparse.cpp
)

add_library(seamass_math ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(seamass_math PUBLIC
  ${MKL_INCLUDE_DIR}
)
target_link_libraries(seamass_math
  ${MKL_LIBRARIES}
)
