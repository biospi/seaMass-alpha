find_package(Boost 1.54 REQUIRED COMPONENTS
	program_options
)

add_executable(seamass
	seamass.cpp
)
target_include_directories(seamass PUBLIC
	${Boost_INCLUDE_DIRS}
)
target_link_libraries (seamass LINK_PUBLIC
	seamass_io
	seamass_core	
	${IntelOpenMP_LIBRARIES}
	${Boost_LIBRARIES}
)
