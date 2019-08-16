find_package(PugiXML REQUIRED)
find_package(TIFF REQUIRED)

message(TIFF_FOUND: ${TIFF_FOUND})
message(TIFF Library: ${TIFF_LIBRARY})
message(TIFF Include: ${TIFF_INCLUDE_DIR})

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
	)
target_link_libraries(seamass_core
        seamass_asrl
        ${PugiXML_LIBRARIES}
		${TIFF_LIBRARY}
	)
