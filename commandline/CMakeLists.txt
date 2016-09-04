find_package(Boost 1.54 COMPONENTS program_options filesystem system REQUIRED)

add_executable(mzmlb2smi
	mzmlb2smi.cpp
)
target_include_directories(mzmlb2smi PUBLIC
	${Boost_INCLUDE_DIRS}
)
target_link_libraries (mzmlb2smi LINK_PUBLIC
	seamass_io
	seamass_core	
	${IntelOpenMP_LIBRARIES}
	${Boost_LIBRARIES}
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

add_executable(seamass2controlpoints
	seamass2controlpoints.cpp
)
target_include_directories(seamass2controlpoints PUBLIC
	${Boost_INCLUDE_DIRS}
)
target_link_libraries (seamass2controlpoints LINK_PUBLIC
	seamass_io
	seamass_core	
	${IntelOpenMP_LIBRARIES}
	${Boost_LIBRARIES}
)
