find_package(Boost 1.54 COMPONENTS program_options REQUIRED)

set(HEADER_FILES
	#core.hpp
	imagecore.hpp
	#SMGWriter.hpp
	SmrWriter.hpp
)

set(SOURCE_FILES
	#core.cpp
	imagecore.cpp
	#SMGWriter.cpp
	smimage.cpp
	SmrWriter.cpp
)

add_executable(smimage ${HEADER_FILES} ${SOURCE_FILES}
)

target_include_directories(seamass PUBLIC
	${Boost_INCLUDE_DIRS}
)
target_link_libraries(smimage LINK_PUBLIC
	seamass_kernel
	${Boost_LIBRARIES}
)
