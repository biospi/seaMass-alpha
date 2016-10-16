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
	Optimizer.hpp
	OptimizerSrl.hpp
	OptimizerAccelerationEve1.hpp
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
	Optimizer.cpp 
	OptimizerSrl.cpp 
	OptimizerAccelerationEve1.cpp
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
