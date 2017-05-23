find_package(PugiXML REQUIRED)

add_library(seamass_core
        Bspline.hpp
        Bspline.cpp
        BasisBspline.hpp
        BasisBspline.cpp
        BasisBsplineMz.hpp
        BasisBsplineMz.cpp
        BasisBsplineScantime.cpp
        BasisBsplineScantime.hpp
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
        )
target_include_directories(seamass_core PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PugiXML_INCLUDE_DIR}
        )
target_link_libraries(seamass_core
        seamass_asrl
        ${PugiXML_LIBRARIES}
        )
