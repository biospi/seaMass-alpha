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
        seamass_io
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
        ${Boost_LIBRARIES}
        )

add_executable(centroid
        centroid.cpp
        )
target_include_directories(centroid PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (centroid LINK_PUBLIC
        seamass_io
        ${Boost_LIBRARIES}
        )

add_executable(smb2mzmlb
        smb2mzmlb.cpp
        )
target_include_directories(smb2mzmlb PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (smb2mzmlb LINK_PUBLIC
        seamass_io
        ${Boost_LIBRARIES}
        )

add_executable(seamass-td
        seamass-td.cpp
        )
target_include_directories(seamass-td PUBLIC
        ${Boost_INCLUDE_DIRS}
        )
target_link_libraries (seamass-td LINK_PUBLIC
        seamass_topdown
        ${Boost_LIBRARIES}
        )