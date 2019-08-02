find_package(Boost 1.54 COMPONENTS program_options filesystem system REQUIRED)
find_package(TIFF REQUIRED)

add_executable(asrl
        asrl.cpp
        )
target_include_directories(asrl PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (asrl LINK_PUBLIC
        seamass_asrl
        ${Boost_LIBRARIES}
        )

add_executable(genisodists
        genisodists.cpp
        )
target_include_directories(genisodists PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (genisodists LINK_PUBLIC
        seamass_core
        ${Boost_LIBRARIES}
        )

add_executable(msp2sml
        msp2sml.cpp
		)
target_include_directories(msp2sml PUBLIC
		${Boost_INCLUDE_DIRS}
		)
target_link_libraries (msp2sml LINK_PUBLIC
		seamass_core
		${Boost_LIBRARIES}
		)

add_executable(mgf2sml
	mgf2sml.cpp
	)
target_include_directories(mgf2sml PUBLIC
	${Boost_INCLUDE_DIRS}
	)
target_link_libraries (mgf2sml LINK_PUBLIC
	seamass_core
	${Boost_LIBRARIES}
	)

add_executable(tiff2sml
		tiff2sml.cpp
		)
target_include_directories(tiff2sml PUBLIC
		${Boost_INCLUDE_DIRS}
		${TIFF_INCLUDE_DIR}
		)
target_link_libraries (tiff2sml LINK_PUBLIC
		seamass_io
		${Boost_LIBRARIES}
		${TIFF_LIBRARY}
		)

add_executable(seamass
        seamass.cpp
        )
target_include_directories(seamass PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (seamass LINK_PUBLIC
        seamass_core
        ${Boost_LIBRARIES}
        )

add_executable(seamass-restore
        seamass-restore.hpp
        seamass-restore.cpp
        )
target_include_directories(seamass-restore PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (seamass-restore LINK_PUBLIC
        seamass_core
        ${Boost_LIBRARIES}
        )

install(TARGETS seamass seamass-restore 
        RUNTIME DESTINATION bin
        )
