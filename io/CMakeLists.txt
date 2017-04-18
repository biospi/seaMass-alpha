find_package(Boost 1.54 COMPONENTS filesystem REQUIRED)
find_package(PugiXML REQUIRED)

add_library(seamass_io
        Dataset.cpp
        DatasetMzmlb.hpp
        mzMLxml.hpp
        DatasetSeamass.hpp
        Dataset.cpp
        DatasetMzmlb.cpp
        mzMLxml.cpp
        DatasetSeamass.cpp
        mzMLxml.tpp
        DatasetMzmlb.tpp
		)
target_include_directories(seamass_io PUBLIC
        ${Boost_INCLUDE_DIR}
        ${PugiXML_INCLUDE_DIR}
        )
target_link_libraries(seamass_io
        seamass_core
        seamass_kernel
        ${Boost_LIBRARY}
        ${PugiXML_LIBRARIES}
        )

