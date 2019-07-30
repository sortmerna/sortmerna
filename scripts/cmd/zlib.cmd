@echo off
::============================================================================
:: FILE: zlib.cmd
:: Created: Jul 29, 2019 Mon
::
:: zlib [release | debug] [t1]
::============================================================================

setlocal

if "%1" == "" ( set btype=Release ) else (
    if "%1" == "release" ( set btype=Release ) else ( set btype=Debug )
)

if "%2" == "" ( set ptype=t1 ) else ( set ptype=%2 )

@echo on
:build
set ZLIB_HOME=%USERPROFILE:\=/%/a03_libs/zlib
:: trim whitespace, otherwise a space will be put between INSTALL_PREFIX and '/bin'. Strange bug.
set CMAKE_INSTALL_PREFIX=%ZLIB_HOME: =%/dist/%ptype: =%
set INSTALL_BIN_DIR=%CMAKE_INSTALL_PREFIX%/bin
set INSTALL_INC_DIR=%CMAKE_INSTALL_PREFIX%/include
set INSTALL_LIB_DIR=%CMAKE_INSTALL_PREFIX%/lib
set INSTALL_MAN_DIR=%CMAKE_INSTALL_PREFIX%/share/man
set INSTALL_PKGCONFIG_DIR=%CMAKE_INSTALL_PREFIX%/share/pkgconfig
set CMAKE_GEN=Visual Studio 16 2019
::set CMAKE_GEN=Visual Studio 15 2017 Win64

:: substitute fwd slash with bkwd
set wpath=%ZLIB_HOME:/=\%
mkdir "%wpath%\build"
pushd "%wpath%\build"

:: print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
::"%VS_HOME%"\bin\Hostx86\x86\cl.exe

cmake -G "%CMAKE_GEN%" ^
-DINSTALL_BIN_DIR=%INSTALL_BIN_DIR% ^
-DINSTALL_INC_DIR=%INSTALL_INC_DIR% ^
-DINSTALL_LIB_DIR=%INSTALL_LIB_DIR% ^
-DINSTALL_MAN_DIR=%INSTALL_MAN_DIR% ^
-DINSTALL_PKGCONFIG_DIR=%INSTALL_PKGCONFIG_DIR% ^
-DCMAKE_INSTALL_PREFIX=%CMAKE_INSTALL_PREFIX%  ..

cmake --build . --config %btype% --target install

popd

endlocal