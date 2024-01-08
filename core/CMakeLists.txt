find_package(Boost 1.54 COMPONENTS filesystem REQUIRED)
find_package(PugiXML REQUIRED)
find_package(TIFF REQUIRED)

add_library(seamass_core
        Bspline.hpp
        Bspline.cpp
        BasisBspline.hpp
        BasisBspline.cpp
        BasisBsplineMz.hpp
        BasisBsplineMz.cpp
        BasisBsplinePeak.hpp
        BasisBsplinePeak.cpp
        BasisBsplineScantime.cpp
        BasisBsplineScantime.hpp
        BasisBsplineCharge.cpp
        BasisBsplineCharge.hpp
        BasisBsplineScale.cpp
        BasisBsplineScale.hpp
        Seamass.cpp
        Seamass.hpp
        Dataset.cpp
        DatasetMzmlb.hpp
        DatasetSeamass.hpp
        Dataset.cpp
        DatasetMzmlb.cpp
        DatasetSeamass.cpp
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
