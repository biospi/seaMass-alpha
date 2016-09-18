set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MKL REQUIRED)

set(HEADER_FILES
	Matrix.hpp
	MatrixSparse.hpp
	BSpline.hpp
	Basis.hpp
	BasisBspline.hpp
	BasisBsplineMz.hpp
	BasisBsplineScantime.hpp
	BasisBsplineScale.hpp
	OptimizerAsrl.hpp
	SeaMass.hpp
)
set(SOURCE_FILES
	Matrix.cpp
	MatrixSparse.cpp
	Bspline.cpp
	Basis.cpp
	BasisBSpline.cpp
	BasisBSplineMz.cpp
	BasisBSplineScantime.cpp
	BasisBSplineScale.cpp
	OptimizerAsrl.cpp
	SeaMass.cpp
)

add_library(seamass_core ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(seamass_core PUBLIC
  ${MKL_INCLUDE_DIR}
)
target_link_libraries(seamass_core
  ${IntelOpenMP_LIBRARIES}
  ${MKL_LIBRARIES}
)
