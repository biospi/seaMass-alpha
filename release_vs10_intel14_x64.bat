@pushd "%~dp0"
@mkdir release
@cd release
@cmake -G"Visual Studio 10 2010 Win64" -T"Intel C++ Compiler XE 14.0" -C"..\..\seamass-windeps\release.cmake" %* ..
@popd