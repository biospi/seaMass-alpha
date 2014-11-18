@set HDF5_DIR=C:\Program Files\HDF_Group\HDF5\1.8.12\cmake\hdf5
@set NETCDF_DIR=C:\Program Files (x86)\netCDF 4.3.2
@mkdir build
@cd build
@cmake -G"Visual Studio 10 2010 Win64" -T"Intel C++ Compiler XE 14.0" ..
@cd ..