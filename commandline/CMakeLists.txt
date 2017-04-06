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
        ${Boost_LIBRARIES}
)

add_executable(seamass-asrl
        seamass-asrl.cpp
)
target_include_directories(seamass-asrl PUBLIC
        ${Boost_INCLUDE_DIRS}
)
target_link_libraries (seamass-asrl LINK_PUBLIC
        seamass_asrl
        ${Boost_LIBRARIES}
)

add_executable(seamass-core
        seamass-core.cpp
)
target_include_directories(seamass-core PUBLIC
        ${Boost_INCLUDE_DIRS}
)
target_link_libraries (seamass-core LINK_PUBLIC
        seamass_io
        seamass_core
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
        ${Boost_LIBRARIES}
)

add_executable(seamass-td
        seamass-td.cpp
)
target_include_directories(seamass-td PUBLIC
        ${Boost_INCLUDE_DIRS}
)
target_link_libraries (seamass-td LINK_PUBLIC
        seamass_io
        seamass_core
        seamass_topdown
        ${Boost_LIBRARIES}
)