# Building from source code

# Contents

* [General notes](#general-notes)
* [Obtaining the source code](#obtaining-the-source-code)
	* [Download GitHub Release](#download-github-release)
	* [Clone the GitHub Repository and use the head or a release](#clone-the-github-repository-and-use-the-head-or-a-release)
* [building on Linux](#building-on-linux)
* [Building on Windows](#building-on-windows)

## General notes

General build steps:
1. Prepare the build environment
2. Get the sources from GitHub or from GitHub releases
3. Build

`HINT`: The minute details on setting the build environment and performing the build can be found in [travis.yml](https://github.com/biocore/sortmerna/blob/master/.travis.yml), which is a [Travis CI](https://travis-ci.org) configuration file used for Sortmerna automated builds.

SortmeRNA 4.0 is C++17 compliant, and requires a compiler that supports `filesystem` standard library.
Currently only GCC-9 and Windows SDK compiler meet that requirement.

Ready to use GCC-9 distribution only available on Debian (Ubuntu). It is not available on Centos, i.e. no Devtoolset-9 yet, and building GCC-9 on Centos is not trivial.

LLVM 9.0 due to be released later this year will have support for `filesystem`. At that time we'll get back to supporting builds with Clang on Linux and OSX.

## Building on Linux

### Install GCC 9

This is for Debian distros (Ubuntu)

```
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install gcc-9 g++-9
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-9
sudo update-alternatives --install /usr/bin/cpp cpp-bin /usr/bin/cpp-9 60
gcc --version
  gcc (Ubuntu 9.1.0-2ubuntu2~18.04) 9.1.0
```

### Install Devtoolset-9

TODO. This section is for RHEL (Centos) distros. Waiting for Devtoolset-9 to become available.

### Install Conda

### get Sortmerna sources
```
# clone the repository
git clone https://github.com/biocore/sortmerna.git

# alternatively get the release sources
wget https://github.com/biocore/sortmerna/archive/v4.0.0.tar.gz

tar xzf v4.0.0.tar.gz

pushd sortmerna

# If you need a particular release (tag)
git checkout v4.0.0
```

### build

Use python script provided with the Sortmerna distro.
By default the build is performed in the User directory.
All required third party libraries will be checked and installed automatically (in User directory by default)
The default build won't interfere with any existing system installation.

By default the build produces statically linked executable i.e. portable.

```
# navigate to the Sortmerna source directory SMR_HOME
SMR_HOME=$HOME/sortmerna
pushd $SMR_HOME

# modify configuration (optional)
vi build/env.yaml

# run the build
python scripts/build.py --name all [--env $SMR_HOME/script/my_env.yaml]
```

## Building on Windows

MS Visual Studio Community edition and CMake for Windows are required for building SortMeRNA.

(1) Download and Install VS Community edition from [Visual Studio community website](https://www.visualstudio.com/vs/community/)

(2) Install CMake

```
python scripts/build.py --name cmake [--env $SMR_HOME/script/my_env.yaml]
set PATH=C:\libs\cmake-3.13.1-win64-x64\bin;%PATH%
cmake --version
cmake-gui
```

(3) Install Git for Windows

Download binary distribution either portable or the installer from [here](https://git-scm.com/download/win)

The portable distribution is a self-extracting archive that can be installed in a directory of your choice e.g.

```
C:\libs\git-2.16.2-64\
	bin\
	cmd\
	dev\
	etc\
	mingw64\
	tmp\
	usr\
```
You can use either `bash.exe` provided with the Git distro, or native Windows CMD `cmd.exe`.

If you choose to work with `cmd.exe`, add the following to your path:

```
set GIT_HOME=C:\libs\git-2.16.2-64
set PATH=%GIT_HOME%\bin;%GIT_HOME%\usr\bin;%GIT_HOME%\mingw64\bin;%PATH%

git --version
```

(4) Configure and build Zlib library

Because we use CMake as the build tool, the steps are the same a for Linux (see Linux instructions above).
For the toolchain specify `cmake -G "Visual Studio 15 2017 Win64"`.
No need to specify `-DCMAKE_BUILD_TYPE` because Visual studio is a multi-type configuration in CMake parlant.

```
git clone https://github.com/madler/zlib.git
```

If using command line the steps are the same as on Linux (see above)

If using CMake GUI:
- click `Browse Source...` and select `%ZLIB_HOME%` i.e. the root directory where the ZLib sources are.
- click `Browse Build...` and select `%ZLIB_HOME%\build\` (confirm to create the `build` directory if not already exists)
- click `Configure` and set the required variables or accept defaults
- click `Generate`

In Visual Studio
- `File -> Open -> Project/Solution` and select `%ZLIB_HOME%\build\zlib.sln`
- In Solution Explorer right-click `ALL_BUILD` and select `build` from drop-down menu

(5) Configure and build RockDB library

```
# modify thirdparty.inc to point to a correct ZLIB installation.
# For example:
...
set(ZLIB_HOME C:/libs/zlib/dist)
set(ZLIB_INCLUDE ${ZLIB_HOME}/include)
set(ZLIB_LIB_DEBUG ${ZLIB_HOME}/lib/zlibstaticd.lib)
set(ZLIB_LIB_RELEASE ${ZLIB_HOME}/lib/zlibstatic.lib)
...
```

If using command line:

```
set SMR_HOME=C:/myprojects/sortmerna
mkdir %SMR_HOME%\build
pushd %SMR_HOME%\build

cmake -G "Visual Studio 15 2017 Win64" \
	-DCMAKE_INSTALL_PREFIX=%SMR_HOME%/dist \
	-DCPACK_BINARY_7Z=ON \
	-DCPACK_SOURCE_7Z=ON \
	-DCPACK_SOURCE_ZIP=OFF \
	-DWITH_MD_LIBRARY=ON \
	-DZLIB_ROOT="C:/libs/zlib/dist" \
	-DZLIB_LIBRARY_RELEASE="C:/libs/zlib/dist/lib/zlibstatic.lib" \
	-DZLIB_LIBRARY_DEBUG="C:/libs/zlib/dist/lib/zlibstaticd.lib" \
	-DROCKSDB_SRC="C:/libs/rocksdb" \
	-DROCKSDB_HOME="C:/libs/rocksdb/dist/d4" \
	-DRAPIDJSON_HOME="C:/libs/rapidjson/dist" \
	-DDIRENTWIN_HOME=C:/libs/dirent ..

cmake --build .                   # build DEBUG (default)
# 		OR
cmake --build . --config Release  # build RELEASE

cmake --build . --target install
cmake --build . --target package
```

If using CMake GUI:
- click `Browse Source...` and select `%ROCKSDB_SRC%"`
- click `Browse Build...` and select `%ROCKSDB_SRC%\build\` (confirm to create the `build` directory if not already exists)
- click `Configure` and set the following variables:
  - Ungrouped Entries
    - PORTABLE (check)
	- GIT_EXECUTABLE (select path to `git.exe` e.g. `C:/libs/git-2.16.2-64/bin/git.exe`
  - WITH
    - WITH_MD_LIBRARY=1
	- WITH_ZLIB=1
	- WITH_TESTS=0
	- WITH_TOOLS=0
	- Accept defaults for the rest
  - CMAKE
	- CMAKE_INSTALL_PREFIX=../../dist
- click `Generate`

In Visual Studio
- `File -> Open -> Project/Solution` and select `%ROCKSDB_SRC%\build\rocksdb.sln`
- In Solution Explorer right-click `ALL_BUILD` and select `build` from drop-down menu

(6) Configure and build RapidJson

See instructions for Linux if using command line.

(7) Configure and build Dirent for Windows [2]

This is optional. If CMake won't find the Dirent when configuring Sortemrna, it will download it from GitHub [2] into `sortmerna/3rdparty/` folder. Only a single header file is used, so no build is required.

(8) Configure and build Sortmerna

Using standard `cmd` Windows shell:

`SMR_HOME` is the top directory where SortMeRNA source distribution (e.g. Git repo) is installed.

```
git clone https://github.com/biocore/sortmerna.git

set SMR_HOME=C:/projects/sortmerna 

mkdir %SMR_HOME%/build
pushd %SMR_HOME%/build

cmake -G "Visual Studio 15 2017 Win64" \
	-DCMAKE_INSTALL_PREFIX=%SMR_HOME%/dist \
	-DCPACK_BINARY_7Z=ON \
	-DCPACK_SOURCE_7Z=ON \
	-DCPACK_SOURCE_ZIP=OFF \
	-DWITH_MD_LIBRARY=ON \
	-DZLIB_ROOT="%LIBS_HOME%/zlib/dist" \
	-DZLIB_LIBRARY_RELEASE="%ZLIB_HOME%/dist/lib/zlibstatic.lib" \
	-DZLIB_LIBRARY_DEBUG="%ZLIB_HOME%/dist/lib/zlibstaticd.lib" \
	-DROCKSDB_SRC="%ROCKSDB_SRC%" \
	-DROCKSDB_HOME="%ROCKSDB_SRC%/dist" \
	-DRAPIDJSON_HOME="%RAPIDJSON_SRC/dist" \
	-DDIRENTWIN_HOME=%DIRENTWIN_SRC%/dist ..

# possible builds: Debug | Release | RelWithDebInfo | MinSizeRel
# Debug - default
cmake --build .
#     OR Release
cmake --build . --config Release

# generate distro
cmake --build . --target install
cmake --build . --target package

# test
set PATH=%SMR_HOME%\dist\bin;%PATH%
sortmerna --version
```

If using CMake GUI

To launch `CMake GUI` navigate to CMake installation directory (using Windows Explorer) and double-click
`cmake-gui`, or launch it from command line:

```
set PATH=C:\libs\cmake-3.13.1-win64-x64\bin;%PATH%
cmake-gui
```

In GUI:
 - click `Browse source` button and navigate to the directory where Sortmerna sources are located (SMR_HOME).
 - click `Browse Build` and navigate to the directory where to build the binaries e.g. `%SMR_HOME%\build`
 - at the prompt select the Generator from the list e.g. "Visual Studio 15 2017 Win64"
 - click `Configure`

 - Set the same variables as shown above for command line i.e.:
	- CMAKE_INSTALL_PREFIX
	- CPACK_BINARY_7Z
	- CPACK_SOURCE_7Z
	- CPACK_SOURCE_ZIP
	- WITH_MD_LIBRARY
	- ZLIB_ROOT
	- ZLIB_LIBRARY_RELEASE
	- ZLIB_LIBRARY_DEBUG
	- ROCKSDB_SRC
	- ROCKSDB_HOME
	- RAPIDJSON_HOME
	- DIRENTWIN_HOME

 - click `Configure` again
 - click `Generate` if all variables are OK (no red background)

The `Generate` generates VS project files in `%SMR_HOME%\build\` directory.
Start Visual Studio and:
- `File -> Open -> Project/Solution .. open %SMR_HOME%\build\sortmerna.sln`
- Select desired build type: `Release | Debug | RelWithDebInfo | MinSizeRel`.
- In Solution explorer right-click `ALL_BUILD' and select `build` in drop-down menu.
- Execute target 'INSTALL'
- Execute target 'PACKAGE'

## Building on Mac

TODO