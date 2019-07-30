@echo off
::============================================================================
:: FILE: rocksdb.cmd
:: Created: Jul 29, 2019 Mon
::
:: build [release | debug] [t3]
::============================================================================

setlocal

if "%1" == "" ( 
    set btype=Release
    set ptype=t3
) else ( if "%2" == "" ( set ptype=t3 ) else ( set ptype=%2 ))
if "%1" == "debug"   set btype=Debug
if "%1" == "release" set btype=Release

if "%ptype%" == "t3" (
    echo "Building Type 3 package: /MD + static ZLib"
)

@echo on
:build
set ROCKSDB_HOME=C:\Users\ak\a03_libs\rocksdb
::set PORTABLE=1
set OPTDBG=0
set WITH_MD_LIBRARY=1
set WITH_RUNTIME_DEBUG=0
set WITH_ZLIB=1
set WITH_XPRESS=1
set WITH_TESTS=0
set WITH_TOOLS=0
set ROCKSDB_INSTALL_ON_WINDOWS=1
set CMAKE_INSTALL_PREFIX=%ROCKSDB_HOME%/dist/%ptype%/%btype%
set CMAKE_GEN=Visual Studio 16 2019
::set CMAKE_GEN=Visual Studio 15 2017 Win64

mkdir %ROCKSDB_HOME%\build
pushd %ROCKSDB_HOME%\build

:: print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
::"%VS_HOME%"\bin\Hostx86\x86\cl.exe

cmake -G "%CMAKE_GEN%" ^
-DOPTDBG=%OPTDBG% ^
-DWITH_MD_LIBRARY=%WITH_MD_LIBRARY% ^
-DWITH_RUNTIME_DEBUG=%WITH_RUNTIME_DEBUG% ^
-DWITH_ZLIB=%WITH_ZLIB% ^
-DWITH_XPRESS=%WITH_XPRESS% ^
-DWITH_TESTS=%WITH_TESTS% ^
-DWITH_TOOLS=%WITH_TOOLS% ^
-ROCKSDB_INSTALL_ON_WINDOWS=%ROCKSDB_INSTALL_ON_WINDOWS% ^
-DCMAKE_INSTALL_PREFIX=%CMAKE_INSTALL_PREFIX%  ..

cmake --build . --config %btype% --target install

popd

endlocal