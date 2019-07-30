@echo off
::============================================================================
:: FILE: rapidjson.cmd
:: Created: Jul 29, 2019 Mon
::
:: build [release | debug]
::============================================================================

setlocal

if "%1" == ""        set btype=Release
if "%1" == "debug"   set btype=Debug
if "%1" == "release" set btype=Release

@echo on
:build
set RAPIDJSON_HOME=C:\Users\ak\a03_libs\rapidjson
set RAPIDJSON_BUILD_EXAMPLES=0
set RAPIDJSON_BUILD_DOC=0
set CMAKE_INSTALL_PREFIX=%RAPIDJSON_HOME%/dist
set CMAKE_GEN=Visual Studio 16 2019
::set CMAKE_GEN=Visual Studio 15 2017 Win64

mkdir %RAPIDJSON_HOME%\build
pushd %RAPIDJSON_HOME%\build

:: print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
::"%VS_HOME%"\bin\Hostx86\x86\cl.exe

cmake -G "%CMAKE_GEN%" ^
-DRAPIDJSON_BUILD_EXAMPLES=%RAPIDJSON_BUILD_EXAMPLES% ^
-DRAPIDJSON_BUILD_DOC=%RAPIDJSON_BUILD_DOC% ^
-DCMAKE_INSTALL_PREFIX=%CMAKE_INSTALL_PREFIX%  ..

cmake --build . --config %btype% --target install

popd

endlocal