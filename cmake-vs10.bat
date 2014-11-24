@set HDF5_DIR=C:\Program Files\HDF_Group\HDF5\1.8.12\cmake\hdf5
@set BOOST_ROOT=..\..\seamass-winlib\boost_1_57_0
@set OSGEO4W_ROOT=..\..\seamass-winlib\spatialindex-src-1.8.5\win64
@mkdir build
@cd build
@cmake -G"Visual Studio 10 2010 Win64" -T"Intel C++ Compiler XE 14.0" -DBoost_USE_STATIC_LIBS=true -DBoost_USE_STATIC=true %* ..
@cd ..