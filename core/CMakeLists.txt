#find_package(OpenMP REQUIRED)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_library(seamass_core seaMassRestoration.cpp OptimiserASRL.cpp HDF5File.cpp BasisFunctions.cpp core.cpp)
target_include_directories(seamass_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS seamass_core DESTINATION lib)
install(FILES seaMassRestoration.hpp DESTINATION include)