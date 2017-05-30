find_package(Boost 1.54 COMPONENTS program_options filesystem system REQUIRED)

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
        seamass-restore.cpp
        )
target_include_directories(seamass-restore PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (seamass-restore LINK_PUBLIC
        seamass_core
        ${Boost_LIBRARIES}
        )

add_executable(seamass-peak
        seamass-peak.cpp
		seamass-peak.hpp
        )
target_include_directories(seamass-peak PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (seamass-peak LINK_PUBLIC
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

install(TARGETS mzmlb2smb seamass seamass-restore seamass-peak smb2mzmlb
        RUNTIME DESTINATION bin
        )
