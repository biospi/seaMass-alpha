@pushd "%~dp0"
@mkdir build
@cd build
@cmake -G"Visual Studio 10 2010 Win64" -T"Intel C++ Compiler XE 14.0" -C"..\..\seamass-windeps\build.cmake" %* ..
@popd.