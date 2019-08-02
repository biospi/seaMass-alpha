find_package(PugiXML REQUIRED)

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
        DatasetSeamass.hpp
        Dataset.cpp
        DatasetSeamass.cpp
        )
target_include_directories(seamass_core PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PugiXML_INCLUDE_DIR}
        )
target_link_libraries(seamass_core
        seamass_asrl
        ${PugiXML_LIBRARIES}
        )
