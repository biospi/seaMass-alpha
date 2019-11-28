find_package(Boost 1.54 COMPONENTS program_options REQUIRED)
find_package(PugiXML REQUIRED)

add_executable(peak
		BsplineData.hpp
		MathOperator.hpp
		PeakData.hpp
		PeakManager.hpp
		PeakOperator.hpp
		SMData.hpp
		smpeak.cpp
		BsplineData.tpp
		MathOperator.tpp
		PeakData.tpp
		PeakOperator.tpp
		SMData.tpp
		)
target_include_directories(peak PUBLIC
		${Boost_INCLUDE_DIRS}
		${PugiXML_INCLUDE_DIR}
		)
target_link_libraries(peak
		seamass_core
		${Boost_LIBRARIES}
		${PugiXML_LIBRARIES}
		)
