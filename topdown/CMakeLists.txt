set(HEADER_FILES
  BasisChargeDistribution.hpp
  BasisIsotopeDistribution.hpp
  SeamassTopDown.hpp
)
set(SOURCE_FILES
  BasisChargeDistribution.cpp
  BasisIsotopeDistribution.cpp
  SeamassTopDown.cpp
)

add_library(seamass_topdown ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(seamass_topdown
  seamass_math
)