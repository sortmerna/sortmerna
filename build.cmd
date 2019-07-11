@echo off
::============================================================================
:: FILE: build.cmd
:: Created: Jul 10, 2019 Wed
::
:: build [release | debug] [t3]
::============================================================================

setlocal

if "%1" == "" ( 
    set btype=Debug
    set ptype=t3
) else ( if "%2" == "" ( set ptype=t3 ) else ( set ptype=%2 ))
if "%1" == "debug"   set btype=Debug
if "%1" == "release" set btype=Release

if "%ptype%" == "t3" (
    echo "Building Type 3 package: /MD + static RocksDB + static ZLib"
)

@echo on
:build
set PORTABLE=0
set WITH_MD_LIBRARY=1
set WITH_RUNTIME_DEBUG=0
set CPACK_BINARY_NSIS=0
set CPACK_BINARY_7Z=1
set CPACK_BINARY_ZIP=1
set CPACK_SOURCE_7Z=1
set CPACK_SOURCE_ZIP=1
set ROCKSDB_STATIC=1
set ZLIB_STATIC=1
set ZLIB_ROOT=C:/a03_libs/zlib/dist/t1
set ZLIB_LIBRARY_RELEASE=%ZLIB_ROOT%/lib/zlibstatic.lib
set ZLIB_LIBRARY_DEBUG=%ZLIB_ROOT%/lib/zlibstaticd.lib
set ROCKSDB_SRC=C:/a03_libs/rocksdb
set ROCKSDB_HOME=C:/a02_projects/rocksdb/dist/Release
set RAPIDJSON_HOME=C:/a03_libs/rapidjson/dist
set DIRENTWIN_HOME=C:/a03_libs/dirent
set CMAKE_INSTALL_PREFIX=C:/a02_projects/sortmerna/dist/%ptype%/%btype%

pushd .\build

:: print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
"%VS_HOME%"\bin\Hostx86\x86\cl.exe

cmake -G "Visual Studio 15 2017 Win64" ^
-DPORTABLE=0 ^
-DCPACK_BINARY_NSIS=0 ^
-DCPACK_BINARY_7Z=1 ^
-DCPACK_BINARY_ZIP=1 ^
-DCPACK_SOURCE_7Z=1 ^
-DCPACK_SOURCE_ZIP=1 ^
-DWITH_MD_LIBRARY=1 ^
-DWITH_RUNTIME_DEBUG=0 ^
-DWITH_TESTS=1 ^
-DZLIB_STATIC=1 ^
-DROCKSDB_STATIC=1 ^
-DCMAKE_INSTALL_PREFIX=%CMAKE_INSTALL_PREFIX% ^
-DZLIB_ROOT=%ZLIB_ROOT% ^
-DZLIB_LIBRARY_RELEASE=%ZLIB_LIBRARY_RELEASE% ^
-DZLIB_LIBRARY_DEBUG=%ZLIB_LIBRARY_DEBUG% ^
-DROCKSDB_SRC=%ROCKSDB_SRC% ^
-DROCKSDB_HOME=%ROCKSDB_HOME% ^
-DRAPIDJSON_HOME=%RAPIDJSON_HOME% ^
-DDIRENTWIN_HOME=%DIRENTWIN_HOME% ..

cmake --build . --config %btype% --target install

%CMAKE_INSTALL_PREFIX%\bin\sortmerna --version
%CMAKE_INSTALL_PREFIX%\bin\sortmerna -h

endlocal