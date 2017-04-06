set(SOURCE_FILES
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
        SeamassCore.cpp
        SeamassCore.hpp
)

add_library(seamass_core ${SOURCE_FILES})
target_link_libraries(seamass_core
    seamass_asrl
)
