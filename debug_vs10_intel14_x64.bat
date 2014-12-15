@pushd "%~dp0"
@mkdir debug
@cd debug
@cmake -G"Visual Studio 10 2010 Win64" -T"Intel C++ Compiler XE 14.0" -C"..\..\seamass-windeps\debug.cmake" %* ..
@popd