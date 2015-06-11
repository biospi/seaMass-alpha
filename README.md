seamass
=======

seaMass sparse signal decomposition and restoration of raw LC-MS data

dependencies
-------
Depends on _CMake_, _HDF5_, _Boost_ and _libspatialindex_. **Note: _libspatialindex_ MUST be at least version 1.8.5**.
Windows developers can get these from the [dependency repository](https://github.com/biospi/seamass-windeps).
Also currently requires _Intel MKL_ with the _Intel C++ Compiler 14.0 or greater_.

building
-------
Run the cmake.sh or cmake.bat file for your platform, which will create make files for compilation in
build/debug and build/release subfolders.
