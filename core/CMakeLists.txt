set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MKL REQUIRED)

set(HEADER_FILES
	BSpline.hpp
	OptimiserASRL.hpp
	BasisFunctions.hpp
	seaMass.hpp
)
set(SOURCE_FILES
	BSpline.cpp
	OptimiserASRL.cpp
	BasisFunctions.cpp
	seaMass.cpp
)

add_library(seamass_core ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(seamass_core PUBLIC
  ${MKL_INCLUDE_DIR}
)
target_link_libraries(seamass_core
  ${IntelOpenMP_LIBRARIES}
  ${MKL_LIBRARIES}
)
