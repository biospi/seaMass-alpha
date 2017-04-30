find_package(Boost 1.54 COMPONENTS program_options REQUIRED)

set(HEADER_FILES
	imagecore.hpp
)

set(SOURCE_FILES
	imagecore.cpp
	smimage.cpp
)

add_executable(smimage ${HEADER_FILES} ${SOURCE_FILES}
)

target_include_directories(smimage PUBLIC
	${Boost_INCLUDE_DIRS}
)
target_link_libraries(smimage LINK_PUBLIC
	seamass_io
	${Boost_LIBRARIES}
)
