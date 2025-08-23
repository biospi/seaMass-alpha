find_package(Boost 1.70 COMPONENTS filesystem REQUIRED NO_MODULE)
find_package(PugiXML REQUIRED)
find_package(TIFF REQUIRED)

add_library(seamass_core
        Bspline.hpp
        Bspline.cpp
        BasisGrid.hpp
        BasisGrid.cpp
        BasisLibrary.cpp
        BasisLibrary.hpp
        BasisBsplineScale.cpp
        BasisBsplineScale.hpp
        BasisBsplinePsf.hpp
        BasisBsplinePsf.cpp
        BasisBsplineScantime.cpp
        BasisBsplineScantime.hpp
        Seamass.cpp
        Seamass.hpp
        Dataset.hpp
        DatasetMzmlb.hpp
        DatasetSeamass.hpp
        DatasetLibrary.hpp
        Dataset.cpp
        DatasetMzmlb.cpp
        DatasetSeamass.cpp
        DatasetLibrary.cpp
        DatasetMzmlb.tpp
        DatasetTiff.cpp
        )
target_include_directories(seamass_core PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PugiXML_INCLUDE_DIR}
        ${TIFF_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
	)
target_link_libraries(seamass_core
        seamass_asrl
        ${PugiXML_LIBRARIES}
        ${TIFF_LIBRARY}
        ${Boost_LIBRARIES}
	)
