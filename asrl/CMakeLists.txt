add_library(seamass_asrl
        Basis.hpp
        Basis.cpp
        BasisMatrix.hpp
        BasisMatrix.cpp
        Optimizer.hpp
        Optimizer.cpp
        OptimizerSrl.hpp
        OptimizerSrl.cpp
        OptimizerAccelerationEve1.hpp
        OptimizerAccelerationEve1.cpp
        SeamassAsrl.hpp
        SeamassAsrl.cpp
        )

target_link_libraries(seamass_asrl
        seamass_kernel
        )
