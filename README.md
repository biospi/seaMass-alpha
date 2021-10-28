seamass
=======

seaMass sparse signal decomposition and restoration of raw Imaging Mass Cytometry (IMC) data.
‘seaMass-IMC’ is a Poisson based sparse spatial regression and Bayesian models for the analysis of large-scale image-based datasets produced by Imaging Mass Cytometry instruments. We have designed seaMass-IMC to process image based ``TIFF`` datasets produced by current IMC instruments and export of denoised TIFF images for downstream processing in IMCs pipelines.
Below is a demonstration of an IMC TIFF image of a biological dataset that has been processed using seaMass-IMC. 

Imaging Mass Cytometry     |  SeaMass Imaging Mass Cytometry
:-------------------------:|:-------------------------:
![CyTof](https://github.com/biospi/seaMass-alpha/blob/feature/libraryImage/plots/142Nd_MHCcII.original.png?raw=true)  |  ![CyTof](https://github.com/biospi/seaMass-alpha/blob/feature/libraryImage/plots/142Nd_MHCcII.seamass2.png?raw=true)


dependencies
-------
Depends on _CMake_, _HDF5_, _Boost_, _libspatialindex_, _pugixml_ and _netCDF4_. **Note: _libspatialindex_ MUST be at least version 1.8.5**.
Windows developers can get these from the [dependency repository](https://github.com/biospi/seamass-windeps).
Also currently requires _Intel MKL_ with the _Intel C++ Compiler 14.0 or 15.0_.

building
-------
Run the cmake.sh or cmake.bat file for your platform, which will create make files for compilation in
build/debug and build/release subfolders.
