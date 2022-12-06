find_package(Boost 1.54 COMPONENTS program_options filesystem system REQUIRED)
find_package(TIFF REQUIRED)

message(TIFF_FOUND: ${TIFF_FOUND})
message(TIFF Library: ${TIFF_LIBRARY})
message(TIFF Include: ${TIFF_INCLUDE_DIR})

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

add_executable(driams_raw2binned
	driams_raw2binned.cpp
	)
target_include_directories(driams_raw2binned PUBLIC
	${Boost_INCLUDE_DIRS}
	)
target_link_libraries (driams_raw2binned LINK_PUBLIC
	seamass_core
	${Boost_LIBRARIES}
	)

add_executable(image2sml
	image2sml.cpp
	)
target_include_directories(image2sml PUBLIC
	${Boost_INCLUDE_DIRS}
	${TIFF_INCLUDE_DIR}
	)
target_link_libraries (image2sml LINK_PUBLIC
	seamass_io
	${Boost_LIBRARIES}
	${TIFF_LIBRARY}
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

add_executable(mzmlb2smb
        mzmlb2smb.cpp
        )
target_include_directories(mzmlb2smb PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (mzmlb2smb LINK_PUBLIC
        seamass_core
        ${Boost_LIBRARIES}
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

add_executable(smb2mzmlb
        smb2mzmlb.cpp
        )
target_include_directories(smb2mzmlb PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries(smb2mzmlb LINK_PUBLIC
        seamass_core
        ${Boost_LIBRARIES}
        )

install(TARGETS mzmlb2smb seamass seamass-restore smb2mzmlb
        RUNTIME DESTINATION bin
        )
